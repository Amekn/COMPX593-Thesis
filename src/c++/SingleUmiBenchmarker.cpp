/*
Usage:
    Benchmarker <ground_truth_fastq> <test1_fastq> <test2_fastq> ...

Optional:
    --band N   (default: 64)  Half-bandwidth for banded global alignment

Description:
    Compares multiple FASTQ files containing processed nanopore reads against a ground truth illumina
    reads FASTQ file. Outputs precision, recall, and F1-score for each test FASTQ file.

    Each read in the ground_truth_fastq must contain single UMI tag in the header
    Each read in the test FASTQ file must contain single UMI tag in the header

    A unordered_multimap is created which allow insertion of multiple pairs of (key, value) pairs with same key.

    Evaluation logic:
        For each test read (UMI key possibly duplicated)
            If UMI key is found in ground truth:
                If key maps to a single reference -> align once
                If key maps to multiple references -> align to all, pick best score, record stats

Important:
    If a ground-truth quality character is '#', that reference position is IGNORED in match/sub/del totals.
    Insertions always count (they don't consume reference positions).

Output (TSV) columns for each test FASTQ:
    Test FASTQ filename
    Number of processed reads (reads compared: those with UMI found in ground truth)
    Number of nucleotides compared (matches + subs + ins + del, with match/sub/del ignoring ref '#')
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

Precision/Recall (nucleotide-level, alignment-based):
    TP = matched
    FP = substitutions + insertions
    FN = substitutions + deletions
    precision = matched / (matched + subs + ins)
    recall    = matched / (matched + subs + del)
*/

#include <algorithm>
#include <cctype>
#include <cstdint>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <vector>

using namespace std;

struct FastqRecord {
    string header;
    string seq;
    string plus;
    string qual;
};

static bool read_fastq_record(istream& in, FastqRecord& r) {
    r = {};
    if (!std::getline(in, r.header)) return false;
    if (!std::getline(in, r.seq))  throw runtime_error("Malformed FASTQ: missing sequence line");
    if (!std::getline(in, r.plus)) throw runtime_error("Malformed FASTQ: missing '+' line");
    if (!std::getline(in, r.qual)) throw runtime_error("Malformed FASTQ: missing quality line");

    if (r.header.empty() || r.header[0] != '@')
        throw runtime_error("Malformed FASTQ: header does not start with '@': " + r.header);
    if (r.plus.empty() || r.plus[0] != '+')
        throw runtime_error("Malformed FASTQ: plus line does not start with '+': " + r.plus);
    if (r.seq.size() != r.qual.size())
        throw runtime_error("Malformed FASTQ: seq/qual length mismatch in record: " + r.header);

    return true;
}

static inline bool is_space(char c) {
    return c == ' ' || c == '\t' || c == '\r' || c == '\n';
}

// Extract UMI after "UMI_" up to whitespace/end.
// Returns empty string if not found.
static string extract_umi_key(const string& header) {
    // header includes leading '@'
    // We search both "UMI_" and "umi_" just in case humans got creative.
    size_t pos = header.find("UMI_");
    if (pos == string::npos) pos = header.find("umi_");
    if (pos == string::npos) return "";

    pos += 4; // move past "UMI_"
    size_t end = pos;
    while (end < header.size() && !is_space(header[end])) end++;
    if (end <= pos) return "";
    return header.substr(pos, end - pos);
}

struct RefRecord {
    string seq;
    string qual; // same length as seq
};

struct AlignmentCounts {
    int score = std::numeric_limits<int>::min()/4;
    uint64_t matches = 0;
    uint64_t subs = 0;
    uint64_t ins = 0;
    uint64_t del = 0;
    uint64_t total = 0; // matches+subs+ins+del (with match/sub/del ignoring ref '#')
};

// Scoring:
//   match: +2
//   mismatch: -1
//   gap: -2
//
// IMPORTANT: if refQual[j] == '#', then
//   - diag contributes 0 (no reward/penalty)
//   - deletion (left move) contributes 0
//   - counts for match/sub/del at that ref position are ignored
static AlignmentCounts banded_global_align(
    const string& query,
    const string& ref,
    const string& refQual,
    int band_half_width
) {
    const int MATCH = 2;
    const int MISMATCH = -1;
    const int GAP = -2;

    const int m = (int)query.size();
    const int n = (int)ref.size();
    if ((int)refQual.size() != n)
        throw runtime_error("Internal error: ref qual length != ref length");

    // Make sure band is at least large enough to connect (m,n) from (0,0)
    // Otherwise dp[m][n] may be unreachable.
    band_half_width = max(band_half_width, abs(m - n) + 2);

    // startJ/endJ for each row i
    vector<int> startJ(m + 1), endJ(m + 1);

    for (int i = 0; i <= m; i++) {
        startJ[i] = max(0, i - band_half_width);
        endJ[i]   = min(n, i + band_half_width);
    }

    // Directions:
    // 0 = diag (consume query+ref)
    // 1 = up   (consume query; insertion relative to ref)
    // 2 = left (consume ref; deletion relative to ref)
    // 255 = start/unused
    vector<vector<uint8_t>> dir(m + 1);

    // Score rows (banded)
    const int NEG_INF = std::numeric_limits<int>::min() / 8;

    vector<int> prevScore, curScore;
    int prevStart = 0;

    // i=0 init
    {
        int i = 0;
        int sJ = startJ[i], eJ = endJ[i];
        int len = eJ - sJ + 1;
        curScore.assign(len, NEG_INF);
        dir[i].assign(len, 255);

        // dp[0][0] = 0 if in band
        if (sJ == 0) {
            curScore[0] = 0;
            dir[i][0] = 255;
        }

        // dp[0][j] from left only
        for (int j = max(1, sJ); j <= eJ; j++) {
            int idx = j - sJ;
            int leftIdx = (j - 1) - sJ;
            int cost = (refQual[j - 1] == '#') ? 0 : GAP; // deleting ref base
            if (leftIdx >= 0 && leftIdx < len && curScore[leftIdx] != NEG_INF) {
                curScore[idx] = curScore[leftIdx] + cost;
                dir[i][idx] = 2; // left
            }
        }

        prevScore.swap(curScore);
        prevStart = sJ;
        curScore.clear();
    }

    // DP for i=1..m
    for (int i = 1; i <= m; i++) {
        int sJ = startJ[i], eJ = endJ[i];
        int len = eJ - sJ + 1;
        curScore.assign(len, NEG_INF);
        dir[i].assign(len, 255);

        // for each j in band
        for (int j = sJ; j <= eJ; j++) {
            int best = NEG_INF;
            uint8_t bestDir = 255;

            // up: dp[i-1][j] + GAP  (insertion relative to ref)
            {
                int pj = j;
                int pIdx = pj - prevStart;
                if (pj >= startJ[i - 1] && pj <= endJ[i - 1]) {
                    int val = prevScore[pIdx];
                    if (val != NEG_INF) {
                        int cand = val + GAP;
                        if (cand > best) {
                            best = cand;
                            bestDir = 1;
                        }
                    }
                }
            }

            // left: dp[i][j-1] + (GAP or 0 if ref '#') (deletion relative to ref)
            if (j > 0) {
                int cj = j - 1;
                if (cj >= sJ) {
                    int cIdx = cj - sJ;
                    int val = curScore[cIdx];
                    if (val != NEG_INF) {
                        int cost = (refQual[j - 1] == '#') ? 0 : GAP;
                        int cand = val + cost;
                        if (cand > best) {
                            best = cand;
                            bestDir = 2;
                        }
                    }
                }
            }

            // diag: dp[i-1][j-1] + match/mismatch/0(if ref '#')
            if (i > 0 && j > 0) {
                int pi = i - 1;
                int pj = j - 1;
                if (pj >= startJ[pi] && pj <= endJ[pi]) {
                    int pIdx = pj - prevStart;
                    int val = prevScore[pIdx];
                    if (val != NEG_INF) {
                        int add = 0;
                        if (refQual[j - 1] != '#') {
                            add = (query[i - 1] == ref[j - 1]) ? MATCH : MISMATCH;
                        }
                        int cand = val + add;
                        // tie-break preference: diag > up > left
                        if (cand > best || (cand == best && bestDir != 0)) {
                            best = cand;
                            bestDir = 0;
                        }
                    }
                }
            }

            int idx = j - sJ;
            curScore[idx] = best;
            dir[i][idx] = bestDir;
        }

        prevScore.swap(curScore);
        prevStart = sJ;
        curScore.clear();
    }

    // If dp[m][n] isn't in band or unreachable, widen band and retry (one retry).
    if (!(n >= startJ[m] && n <= endJ[m])) {
        return banded_global_align(query, ref, refQual, max(m, n));
    }
    int endIdx = n - startJ[m];
    int finalScore = prevScore[endIdx];
    if (finalScore == NEG_INF) {
        return banded_global_align(query, ref, refQual, max(m, n));
    }

    // Backtrace for counts
    AlignmentCounts out;
    out.score = finalScore;

    int i = m, j = n;
    while (i > 0 || j > 0) {
        if (!(j >= startJ[i] && j <= endJ[i])) {
            throw runtime_error("Backtrace stepped خارج band. This should not happen.");
        }
        uint8_t d = dir[i][j - startJ[i]];
        if (d == 0) { // diag
            // consumes query[i-1], ref[j-1]
            char rq = refQual[j - 1];
            if (rq != '#') {
                out.total++;
                if (query[i - 1] == ref[j - 1]) out.matches++;
                else out.subs++;
            }
            i--; j--;
        } else if (d == 1) { // up = insertion relative to ref
            out.total++;
            out.ins++;
            i--;
        } else if (d == 2) { // left = deletion relative to ref
            char rq = refQual[j - 1];
            if (rq != '#') {
                out.total++;
                out.del++;
            }
            j--;
        } else {
            // If we land here, something's wrong with initialization/tie-breaking.
            throw runtime_error("Backtrace hit invalid direction at i=" + to_string(i) + " j=" + to_string(j));
        }
    }

    return out;
}

struct Totals {
    uint64_t reads_compared = 0;
    uint64_t total = 0;
    uint64_t matches = 0;
    uint64_t del = 0;
    uint64_t ins = 0;
    uint64_t subs = 0;
};

static void add_counts(Totals& t, const AlignmentCounts& a) {
    t.total   += a.total;
    t.matches += a.matches;
    t.del     += a.del;
    t.ins     += a.ins;
    t.subs    += a.subs;
}

static double safe_div(double num, double den) {
    return (den == 0.0) ? 0.0 : (num / den);
}

int main(int argc, char** argv) {
    ios::sync_with_stdio(false);
    cin.tie(nullptr);

    int band = 64;

    vector<string> positional;
    for (int i = 1; i < argc; i++) {
        string a = argv[i];
        if (a == "--band") {
            if (i + 1 >= argc) {
                cerr << "Error: --band requires an integer\n";
                return 1;
            }
            band = stoi(argv[++i]);
            if (band < 0) band = 0;
            continue;
        }
        positional.push_back(a);
    }

    if (positional.size() < 2) {
        cerr << "Usage: Benchmarker <ground_truth_fastq> <test1_fastq> <test2_fastq> ... [--band N]\n";
        return 1;
    }

    const string ground_path = positional[0];
    vector<string> test_paths(positional.begin() + 1, positional.end());

    // Load ground truth into unordered_multimap
    unordered_multimap<string, RefRecord> gt;
    gt.reserve(1 << 20);

    {
        ifstream fin(ground_path);
        if (!fin) {
            cerr << "Error: cannot open ground truth FASTQ: " << ground_path << "\n";
            return 1;
        }

        FastqRecord r;
        uint64_t loaded = 0, skipped = 0;
        while (read_fastq_record(fin, r)) {
            string umi = extract_umi_key(r.header);
            if (umi.empty()) {
                skipped++;
                continue;
            }
            gt.emplace(umi, RefRecord{r.seq, r.qual});
            loaded++;
        }

        if (loaded == 0) {
            cerr << "Error: no usable ground truth reads loaded (UMI tag missing?)\n";
            return 1;
        }
        // soft warning, not fatal
        if (skipped > 0) {
            cerr << "[warn] ground truth: skipped " << skipped << " reads with no UMI tag\n";
        }
    }

    // Print TSV header
    cout
        << "Test_FASTQ\t"
        << "Processed_Reads\t"
        << "Nucleotides_Compared\t"
        << "Matched\t"
        << "Deleted\t"
        << "Inserted\t"
        << "Substituted\t"
        << "Norm_Matched\t"
        << "Norm_Deleted\t"
        << "Norm_Inserted\t"
        << "Norm_Substituted\t"
        << "Bases_Mapped_CIGAR\t"
        << "MM2_Mismatch_Rate\t"
        << "Precision\t"
        << "Recall\t"
        << "F1\n";

    cout << fixed << setprecision(6);

    // Process each test file
    for (const auto& test_path : test_paths) {
        ifstream fin(test_path);
        if (!fin) {
            cerr << "Error: cannot open test FASTQ: " << test_path << "\n";
            return 1;
        }

        Totals totals;
        FastqRecord tr;

        uint64_t skipped_no_umi = 0;
        uint64_t skipped_no_gt  = 0;

        while (read_fastq_record(fin, tr)) {
            string umi = extract_umi_key(tr.header);
            if (umi.empty()) {
                skipped_no_umi++;
                continue;
            }

            auto range = gt.equal_range(umi);
            if (range.first == range.second) {
                skipped_no_gt++;
                continue;
            }

            // Align to one or multiple GT candidates and pick best score
            AlignmentCounts best;
            bool have_best = false;

            for (auto it = range.first; it != range.second; ++it) {
                const auto& refRec = it->second;
                AlignmentCounts cur = banded_global_align(tr.seq, refRec.seq, refRec.qual, band);

                if (!have_best || cur.score > best.score) {
                    best = cur;
                    have_best = true;
                }
            }

            if (have_best) {
                totals.reads_compared++;
                add_counts(totals, best);
            }
        }

        double total = (double)totals.total;
        double norm_m = safe_div((double)totals.matches, total);
        double norm_d = safe_div((double)totals.del, total);
        double norm_i = safe_div((double)totals.ins, total);
        double norm_s = safe_div((double)totals.subs, total);

        double bases_mapped_cigar = (double)totals.matches + (double)totals.subs;
        double mm2_mismatch_rate  = safe_div((double)totals.subs, bases_mapped_cigar);

        double precision = safe_div((double)totals.matches,
                                    (double)totals.matches + (double)totals.subs + (double)totals.ins);
        double recall    = safe_div((double)totals.matches,
                                    (double)totals.matches + (double)totals.subs + (double)totals.del);
        double f1 = (precision + recall == 0.0) ? 0.0 : (2.0 * precision * recall / (precision + recall));

        cout
            << test_path << "\t"
            << totals.reads_compared << "\t"
            << totals.total << "\t"
            << totals.matches << "\t"
            << totals.del << "\t"
            << totals.ins << "\t"
            << totals.subs << "\t"
            << norm_m << "\t"
            << norm_d << "\t"
            << norm_i << "\t"
            << norm_s << "\t"
            << (uint64_t)bases_mapped_cigar << "\t"
            << mm2_mismatch_rate << "\t"
            << precision << "\t"
            << recall << "\t"
            << f1 << "\n";

        // Optional stderr warnings (kept out of TSV because humans love breaking parsers)
        if (skipped_no_umi > 0) {
            cerr << "[warn] " << test_path << ": skipped " << skipped_no_umi << " reads with no UMI tag\n";
        }
        if (skipped_no_gt > 0) {
            cerr << "[warn] " << test_path << ": skipped " << skipped_no_gt << " reads with UMI not in ground truth\n";
        }
    }

    return 0;
}
