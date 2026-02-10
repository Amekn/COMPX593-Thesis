/*
Usage: Benchmarker <ground_truth_fastq> <test_fastq1> <test_fastq2> ...

Description:
    Compares multiple FASTQ files containing processed nanopore reads against a ground truth illumina
    reads FASTQ file. Outputs precision, recall, and F1-score for each test FASTQ file.

    Each read in the ground_truth_fastq must contain dual UMI tags in the header
    Each read in the test FASTQ file must contain dual UMI tags in the header (same format as ground truth)
    A map is created from the ground truth FASTQ file using the dual UMI tags as keys.

    Each read in the test FASTQ file is looked up in the map using the dual UMI tags.
    If a match is found, the read sequences are compared to determine statistics of correctness.

Important:
    If a ground-truth quality character is '#', that reference position is IGNORED in match/sub/del totals.

Output (TSV) columns for each test FASTQ:
    Test FASTQ filename
    Number of processed reads
    Number of nucleotides compared
    Number of matched nucleotides
    Number of deleted nucleotides
    Number of inserted nucleotides
    Number of substituted nucleotides
    Normalized Matched Nucleotides (matched / total)
    Normalized Deletions (deletions / total)
    Normalized Insertions (insertions / total)
    Normalized Substitutions (substitutions / total)
    Precision
    Recall
    F1-score

Precision/Recall definition (nucleotide-level, alignment-based):
    TP = matched
    FP = substitutions + insertions
    FN = substitutions + deletions
    precision = TP / (TP + FP) = matched / (matched + subs + ins)
    recall    = TP / (TP + FN) = matched / (matched + subs + del)
*/

#include <algorithm>
#include <cctype>
#include <cstdint>
#include <fstream>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <vector>

using namespace std;

struct GTRead {
    string seq;
    string qual; // same length as seq (after padding/trimming)
};

struct Stats {
    uint64_t processed_reads = 0;
    uint64_t total_compared = 0;   // reference positions considered (qual != '#')
    uint64_t matched = 0;
    uint64_t deleted = 0;          // ref base aligned to gap in test
    uint64_t inserted = 0;         // test base aligned to gap in ref
    uint64_t substituted = 0;      // ref base != test base (both present)
};

static inline char up(char c) {
    return static_cast<char>(toupper(static_cast<unsigned char>(c)));
}

static bool read_fastq_record(istream& in, string& header, string& seq, string& plus, string& qual) {
    header.clear(); seq.clear(); plus.clear(); qual.clear();
    if (!std::getline(in, header)) return false;
    if (!std::getline(in, seq))    throw runtime_error("Malformed FASTQ: missing sequence line");
    if (!std::getline(in, plus))   throw runtime_error("Malformed FASTQ: missing '+' line");
    if (!std::getline(in, qual))   throw runtime_error("Malformed FASTQ: missing quality line");
    if (header.empty() || header[0] != '@')
        throw runtime_error("Malformed FASTQ: header does not start with '@': " + header);
    if (plus.empty() || plus[0] != '+')
        throw runtime_error("Malformed FASTQ: third line does not start with '+': " + plus);
    if (seq.size() != qual.size())
        throw runtime_error("Malformed FASTQ: sequence and quality lengths differ: " + to_string(seq.size()) + " vs " + to_string(qual.size()));
    return true;
}

// Parse ":UMI_<FWD>_<REV>" from header (either with or without the leading '@').
static string parse_dual_umi_key(const string& header_line) {
    const string& h = header_line;
    size_t pos = h.find(":UMI_");
    if (pos == string::npos) {
        throw runtime_error("Header missing ':UMI_' tag: " + h);
    }
    pos += 5; // skip ":UMI_"
    if (pos >= h.size()) {
        throw runtime_error("Header has ':UMI_' but no UMI content: " + h);
    }

    // The UMI segment ends at whitespace (if any)
    size_t end = h.find_first_of(" \t\r\n", pos);
    if (end == string::npos) end = h.size();

    string umi_part = h.substr(pos, end - pos); // "<FWD>_<REV>"
    size_t sep = umi_part.find('_');
    if (sep == string::npos) {
        throw runtime_error("UMI tag not dual (missing underscore): " + h);
    }

    string fwd = umi_part.substr(0, sep);
    string rev = umi_part.substr(sep + 1);
    if (fwd.empty() || rev.empty()) {
        throw runtime_error("Empty UMI component in header: " + h);
    }

    return fwd + "_" + rev; // canonical key
}

// Banded global alignment (Needleman-Wunsch Algorithm)
// with levenshtein costs: match 0, sub 1, ins 1, del 1
// (i.e., editing distance)
// Returns alignment-derived counts; ignores reference positions where ref_qual[i] == '#'.
static void align_and_count(const string& ref, const string& ref_qual, const string& test, Stats& S) {
    // ref: ground truth sequence
    const int n = static_cast<int>(ref.size());
    // test: processed read sequence
    const int m = static_cast<int>(test.size());

    if (n == 0 && m == 0) return;

    const int max_len = max(n, m);
    // Default band: enough to allow typical indels but avoid full O(n*m) for every read.
    int band = max(50, abs(n - m) + 50);
    band = min(band, max_len);

    // Banded DP Lambda Function
    auto run_banded = [&](int use_band) -> bool {
        const uint16_t INF = std::numeric_limits<uint16_t>::max() / 4;

        // There will be (n+1) rows, each with a variable number of columns depending on the band.
        // Each band represents a slice of the test sequence indices valid for that ref index.
        // This assumption is based on the idea that best alignment will not stray too far from the diagonal.
        vector<int> start(n + 1), end(n + 1);
        // dp is the cost matrix; dir is the traceback direction matrix for reconstructing the alignment.
        vector<vector<uint16_t>> dp(n + 1);
        vector<vector<uint8_t>> dir(n + 1); // 0 diag, 1 up, 2 left

        // Row 0 initialization
        start[0] = 0;
        // set band to either default band or length of test seq, whichever is smaller
        end[0] = min(m, use_band);
        // assign a vector of size of the band, initialized to INF
        dp[0].assign(end[0] - start[0] + 1, INF);
        // assign a vector of size of the band, initialized to 255
        dir[0].assign(end[0] - start[0] + 1, 255);
        // initialize first row (i=0)
        for (int j = start[0]; j <= end[0]; ++j) {
            // assign each cell the cost of j
            dp[0][j - start[0]] = static_cast<uint16_t>(j);
            // traceback direction is always from left (insertions)
            // so the first cell (0,0) is a special case
            if (j == 0) dir[0][j - start[0]] = 255;
            // other cells should traceback to the left
            else dir[0][j - start[0]] = 2; // from left (insertions)
        }

        // Fill DP and Direction Matrices
        for (int i = 1; i <= n; ++i) {
            // determine band index range for this row
            start[i] = max(0, i - use_band);
            end[i] = min(m, i + use_band);
            // initialize dp vector for this row using band size
            dp[i].assign(end[i] - start[i] + 1, INF);
            // initialize dir vector for this row using band size
            dir[i].assign(end[i] - start[i] + 1, 255);

            // Fill cells within the band
            for (int j = start[i]; j <= end[i]; ++j) {
                // Compute costs for up (deletion), left (insertion), diag (match/substitution)
                uint16_t best = INF;
                // record best direction: 0 diag, 1 up, 2 left
                uint8_t best_dir = 255;

                // up: (i-1, j) + 1
                if (j >= start[i - 1] && j <= end[i - 1]) {
                    uint16_t up_cost = static_cast<uint16_t>(dp[i - 1][j - start[i - 1]] + 1);
                    best = up_cost;
                    best_dir = 1;
                }

                // left: (i, j-1) + 1
                if (j - 1 >= start[i] && j - 1 <= end[i]) {
                    uint16_t left_cost = static_cast<uint16_t>(dp[i][(j - 1) - start[i]] + 1);
                    if (left_cost < best) {
                        best = left_cost;
                        best_dir = 2;
                    }
                }

                // diag: (i-1, j-1) + (match?0:1)
                if (j - 1 >= 0 && j - 1 >= start[i - 1] && j - 1 <= end[i - 1]) {
                    uint16_t diag_prev = dp[i - 1][(j - 1) - start[i - 1]];
                    uint16_t add = (up(ref[i - 1]) == up(test[j - 1])) ? 0 : 1;
                    uint16_t diag_cost = static_cast<uint16_t>(diag_prev + add);

                    // Tie-break: prefer diag, then up, then left (more stable for substitutions vs indels)
                    if (diag_cost < best || (diag_cost == best && best_dir != 0)) {
                        best = diag_cost;
                        best_dir = 0;
                    }
                }

                dp[i][j - start[i]] = best;
                dir[i][j - start[i]] = best_dir;
            }
        }

        // Ensure target cell (buttom right corner) is inside the band
        if (m < start[n] || m > end[n]) return false;

        // Traceback
        int i = n, j = m;
        while (i > 0 || j > 0) {
            if (j < start[i] || j > end[i]) return false; // outside band during traceback
            uint8_t d = dir[i][j - start[i]];

            if (i == 0) d = 2;         // only insertions
            else if (j == 0) d = 1;    // only deletions

            if (d == 0) { // diag
                char r = up(ref[i - 1]);
                char t = up(test[j - 1]);
                char q = ref_qual[i - 1];
                if (q != '#') {
                    S.total_compared++;
                    if (r == t) S.matched++;
                    else S.substituted++;
                }
                --i; --j;
            } else if (d == 1) { // up: deletion in test (gap in test)
                char q = ref_qual[i - 1];
                if (q != '#') {
                    S.total_compared++;
                    S.deleted++;
                }
                --i;
            } else if (d == 2) { // left: insertion in test
                S.inserted++;
                --j;
            } else {
                // should never happen, but humans love surprises
                return false;
            }
        }

        return true;
    };

    // Try banded; if the band was too tight for some nasty read, fall back to full DP band.
    if (!run_banded(band)) {
        if(!run_banded(max_len)){
            throw runtime_error("Alignment failed even with full DP matrix.");
        }
    }
}

static unordered_map<string, GTRead> load_ground_truth(const string& path) {
    ifstream in(path);
    if (!in) throw runtime_error("Failed to open ground truth FASTQ: " + path);

    unordered_map<string, GTRead> gt;
    // The illumina ground truth containing ~1.7 m reads; reserve to avoid rehashes.
    gt.reserve(3000000);
    gt.max_load_factor(0.70f);

    string h, s, p, q;
    while (read_fastq_record(in, h, s, p, q)) {
        string key = parse_dual_umi_key(h);
        auto it = gt.find(key);
        if (it == gt.end()) {
            gt.emplace(key, GTRead{std::move(s), std::move(q)});
        }
        // If duplicate UMI found in ground truth, ignore later occurrences.
    }
    return gt;
}

static void print_header() {
    cout
        << "test_fastq"
        << "\tprocessed_reads"
        << "\tnucleotides_compared"
        << "\tmatched"
        << "\tdeleted"
        << "\tinserted"
        << "\tsubstituted"
        << "\tnorm_matched"
        << "\tnorm_deleted"
        << "\tnorm_inserted"
        << "\tnorm_substituted"
        << "\tprecision"
        << "\trecall"
        << "\tf1"
        << "\n";
}

static void report(const string& test_path, const Stats& S) {
    double total = (S.total_compared == 0) ? 0.0 : static_cast<double>(S.total_compared);

    double norm_matched = (total == 0.0) ? 0.0 : (double)S.matched / total;
    double norm_del     = (total == 0.0) ? 0.0 : (double)S.deleted / total;
    double norm_ins     = (total == 0.0) ? 0.0 : (double)S.inserted / total;
    double norm_sub     = (total == 0.0) ? 0.0 : (double)S.substituted / total;

    double prec_den = (double)S.matched + (double)S.substituted + (double)S.inserted;
    double rec_den  = (double)S.matched + (double)S.substituted + (double)S.deleted;

    double precision = (prec_den == 0.0) ? 0.0 : (double)S.matched / prec_den;
    double recall    = (rec_den  == 0.0) ? 0.0 : (double)S.matched / rec_den;
    double f1        = ((precision + recall) == 0.0) ? 0.0 : (2.0 * precision * recall) / (precision + recall);

    cout
        << test_path
        << "\t" << S.processed_reads
        << "\t" << S.total_compared
        << "\t" << S.matched
        << "\t" << S.deleted
        << "\t" << S.inserted
        << "\t" << S.substituted
        << "\t" << norm_matched
        << "\t" << norm_del
        << "\t" << norm_ins
        << "\t" << norm_sub
        << "\t" << precision
        << "\t" << recall
        << "\t" << f1
        << "\n";
}

int main(int argc, char** argv) {
    ios::sync_with_stdio(false);
    cin.tie(nullptr);

    try {
        if (argc < 3) {
            cerr << "Usage: " << argv[0] << " <ground_truth_fastq> <test_fastq1> <test_fastq2> ...\n";
            return 1;
        }

        string gt_path = argv[1];
        auto gt = load_ground_truth(gt_path);

        print_header();

        for (int a = 2; a < argc; ++a) {
            string test_path = argv[a];
            ifstream in(test_path);
            if (!in) throw runtime_error("Failed to open test FASTQ: " + test_path);

            Stats S;
            string h, s, p, q;
            uint64_t umi_ok = 0, umi_bad = 0, hits = 0;
            while (read_fastq_record(in, h, s, p, q)) {
                S.processed_reads++;
                string key;
                try {
                    key = parse_dual_umi_key(h);
                    umi_ok++;
                } catch (...) {
                    // Malformed or missing UMI tag in test read: skip this read.
                    umi_bad++;
                    continue;
                }

                auto it = gt.find(key);
                if (it == gt.end()) continue; // no ground-truth match for this UMI

                // Pass ground truth seq/qual and test seq to align_and_count
                // This function updates S in place.
                align_and_count(it->second.seq, it->second.qual, s, S);
            }

            cerr << test_path
                    << " umi_ok=" << umi_ok
                    << " umi_bad=" << umi_bad
                    << " hits=" << hits
                    << " gt_size=" << gt.size()
                    << "\n";
            // Output report for this test FASTQ
            report(test_path, S);
        }

        return 0;
    } catch (const exception& e) {
        cerr << "Error: " << e.what() << "\n";
        return 2;
    }
}
