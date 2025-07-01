#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>
#include <random>
#include <algorithm>
#include <unordered_set>
#include <unordered_map>
#include <cctype>
#include <stdexcept>

using namespace std;

struct Region {
    size_t start;   // inclusive
    size_t end;     // exclusive (must satisfy start < end)
};

// ──────────────────────────────────────────────────────────────────────────────
// Native barcode lookup table (SQK‑NBD114.24, NB01‑NB24)
// Each entry maps the barcode ID (e.g. "NB01") to a pair of
// (forward_barcode, reverse_barcode) sequences.
// Forward strand flanking: 5'‑AAGGTTAA‑barcode‑CAGCACCT‑3'
// Reverse complement strand flanking: 5'‑GGTGCTG‑barcode‑TTAACCTTAGCAAT‑3'
// ──────────────────────────────────────────────────────────────────────────────
static const unordered_map<string, pair<string,string>> BARCODE_SEQS = {
    {"NB01", {"CACAAAGACACCGACAACTTTCTT", "AAGAAAGTTGTCGGTGTCTTTGTG"}},
    {"NB02", {"ACAGACGACTACAAACGGAATCGA", "TCGATTCCGTTTGTAGTCGTCTGT"}},
    {"NB03", {"CCTGGTAACTGGGACACAAGACTC", "GAGTCTTGTGTCCCAGTTACCAGG"}},
    {"NB04", {"TAGGGAAACACGATAGAATCCGAA", "TTCGGATTCTATCGTGTTTCCCTA"}},
    {"NB05", {"AAGGTTACACAAACCCTGGACAAG", "CTTGTCCAGGGTTTGTGTAACCTT"}},
    {"NB06", {"GACTACTTTCTGCCTTTGCGAGAA", "TTCTCGCAAAGGCAGAAAGTAGTC"}},
    {"NB07", {"AAGGATTCATTCCCACGGTAACAC", "GTGTTACCGTGGGAATGAATCCTT"}},
    {"NB08", {"ACGTAACTTGGTTTGTTCCCTGAA", "TTCAGGGAACAAACCAAGTTACGT"}},
    {"NB09", {"AACCAAGACTCGCTGTGCCTAGTT", "AACTAGGCACAGCGAGTCTTGGTT"}},
    {"NB10", {"GAGAGGACAAAGGTTTCAACGCTT", "AAGCGTTGAAACCTTTGTCCTCTC"}},
    {"NB11", {"TCCATTCCCTCCGATAGATGAAAC", "GTTTCATCTATCGGAGGGAATGGA"}},
    {"NB12", {"TCCGATTCTGCTTCTTTCTACCTG", "CAGGTAGAAAGAAGCAGAATCGGA"}},
    {"NB13", {"AGAACGACTTCCATACTCGTGTGA", "TCACACGAGTATGGAAGTCGTTCT"}},
    {"NB14", {"AACGAGTCTCTTGGGACCCATAGA", "TCTATGGGTCCCAAGAGACTCGTT"}},
    {"NB15", {"AGGTCTACCTCGCTAACACCACTG", "CAGTGGTGTTAGCGAGGTAGACCT"}},
    {"NB16", {"CGTCAACTGACAGTGGTTCGTACT", "AGTACGAACCACTGTCAGTTGACG"}},
    {"NB17", {"ACCCTCCAGGAAAGTACCTCTGAT", "ATCAGAGGTACTTTCCTGGAGGGT"}},
    {"NB18", {"CCAAACCCAACAACCTAGATAGGC", "GCCTATCTAGGTTGTTGGGTTTGG"}},
    {"NB19", {"GTTCCTCGTGCAGTGTCAAGAGAT", "ATCTCTTGACACTGCACGAGGAAC"}},
    {"NB20", {"TTGCGTCCTGTTACGAGAACTCAT", "ATGAGTTCTCGTAACAGGACGCAA"}},
    {"NB21", {"GAGCCTCTCATTGTCCGTTCTCTA", "TAGAGAACGGACAATGAGAGGCTC"}},
    {"NB22", {"ACCACTGCCATGTATCAAAGTACG", "CGTACTTTGATACATGGCAGTGGT"}},
    {"NB23", {"CTTACTACCCAGTGAACCTCCTCG", "CGAGGAGGTTCACTGGGTAGTAAG"}},
    {"NB24", {"GCATAGTTCTGCATGATGGGTTAG", "CTAACCCATCATGCAGAACTATGC"}}
};

// Read first sequence from FASTA. Returns sequence; set sequence name
static string read_fasta(const string &path, string &sequence_name) {
    ifstream in(path);
    if (!in)
        throw runtime_error("Cannot open reference FASTA: " + path);
    string line, seq;
    while (getline(in, line)) {
        if (line.empty())
            continue;
        if (line[0] == '>') {
            sequence_name = line.substr(1);
        } else {
            seq += line;
        }
    }
    if (sequence_name.empty() || seq.empty())
        throw runtime_error("Invalid FASTA: " + path);
    return seq;
}

// Read region declaration file (CSV: start,end), return allowed mutable positions
static vector<size_t> read_regions(const string &path, size_t ref_len) {
    ifstream in(path);
    if (!in)
        throw runtime_error("Cannot open region file: " + path);
    vector<size_t> allowed;
    string line;
    while (std::getline(in, line)) {
        if (line.empty())
            continue;
        size_t comma = line.find(',');
        if (comma == string::npos)
            throw runtime_error("Region line missing comma: " + line);
        size_t start = stoul(line.substr(0, comma));
        size_t end   = stoul(line.substr(comma + 1));
        
        if (start >= end)
            throw runtime_error("Region end must be > start: " + line);
        if (end > ref_len)
            throw runtime_error("Region exceeds reference length: " + line);

        for (size_t i = start; i < end; ++i)
            allowed.push_back(i);
    }
    return allowed;
}

static char random_base(char exclude, mt19937 &rng) {
    static const char bases[] = {'A', 'C', 'G', 'T'};
    uniform_int_distribution<int> dist(0, 3);
    char b;
    do {
        b = bases[dist(rng)];
    } while (toupper(b) == toupper(exclude));
    return b;
}

static string random_bases(size_t n, mt19937 &rng) {
    static const char bases[] = {'A', 'C', 'G', 'T'};
    uniform_int_distribution<int> dist(0, 3);
    string s;
    s.reserve(n);
    for (size_t i = 0; i < n; ++i)
        s.push_back(bases[dist(rng)]);
    return s;
}

static inline char comp(char b) {
    switch (toupper(b)) {
        case 'A': return 'T';
        case 'C': return 'G';
        case 'G': return 'C';
        case 'T': return 'A';
        default:  return 'N';
    }
}

static string reverse_complement(const string &seq) {
    string rc;
    rc.reserve(seq.size());
    for (auto it = seq.rbegin(); it != seq.rend(); ++it)
        rc.push_back(comp(*it));
    return rc;
}

// ──────────────────────────────────────────────────────────────────────────────
// main
// ──────────────────────────────────────────────────────────────────────────────

int main(int argc, char *argv[]) {
    // Expected arguments:
    //   1  <variable_mutant_prob>   (0‑100)
    //   2  <framework_mutant_prob>  (0‑100)
    //   3  <max_length_flunctuation>
    //   4  <num_variants>
    //   5  <region_file>
    //   6  <reference_fasta>
    //   7  <output_fasta>
    //   8 [--barcode|-b] (optional)
    //   9 [--revcomp|-r]  (optional)

    bool revcomp_flag  = false;
    bool barcode_flag  = false;   

    // ------------------------------------------------------------------
    // Parse positional arguments (first 7 arguments are mandatory)
    // ------------------------------------------------------------------
    if (argc < 8) {
        cerr << "Usage: " << argv[0]
             << " <variable_mutant_prob> <framework_mutant_prob> <max_length_flunctuation>"
                " <num_variants> <region_file> <reference_fasta> <output_fasta>"
                " [--barcode|-b] [--revcomp|-r]\n";
        return 1;
    }

    // Mandatory positionals
    const size_t var_prob        = stoul(argv[1]);   // 0‑100
    const size_t fw_prob         = stoul(argv[2]);   // 0‑100
    const size_t max_len_flunc   = stoul(argv[3]);
    const size_t num_variants    = stoul(argv[4]);
    const string region_file     = argv[5];
    const string ref_fasta       = argv[6];
    const string out_fasta       = argv[7];

    // Optional flags start at argv[8]
    for (int i = 8; i < argc; ++i) {
        string flag = argv[i];
        if (flag == "--revcomp" || flag == "-r") {
            //Produces variant in reverse complement
            revcomp_flag = true;
        } else if (flag == "--barcode" || flag == "-b") {
            //Produces variant by looping through all barcodes
            barcode_flag = true;         
        } else {
            //If any unknow flag is capture, print an error message.
            cerr << "Unknown flag: " << flag << '\n';
            return 1;
        }
    }

    // Build a sorted list of all available barcodes for easy cycling
    vector<string> all_barcode_ids;
    all_barcode_ids.reserve(BARCODE_SEQS.size());
    for (const auto &kv : BARCODE_SEQS)
        all_barcode_ids.push_back(kv.first);
    sort(all_barcode_ids.begin(), all_barcode_ids.end());   // NB01 .. NB24

    // Check to ensure mutation probability is within the range.
    if (var_prob < 0 || fw_prob < 0 || var_prob > 100 || fw_prob > 100) {
        cerr << "Mutation probabilities must be in the range 0‑100\n";
        return 1;
    }

    // Reference
    string ref_name;
    string ref_seq = read_fasta(ref_fasta, ref_name);
    const size_t ref_len = ref_seq.size();

    // Mutable positions (variable region)
    vector<size_t> allowed = read_regions(region_file, ref_len);
    if (allowed.empty()) {
        cerr << "No regions specified\n";
        return 1;
    }

    // Build mask for quick lookup
    vector<char> is_variable(ref_len, 0);
    for (size_t p : allowed) is_variable[p] = 1;

    // RNG helpers
    mt19937 rng(random_device{}());
    uniform_int_distribution<int> prob_dist(1, 100);  // inclusive
    uniform_int_distribution<int> len_flunc_dist(-static_cast<int>(max_len_flunc), static_cast<int>(max_len_flunc));

    ofstream out(out_fasta);
    if (!out) {
        cerr << "Cannot open output file: " << out_fasta << "\n";
        return 1;
    }

    for (size_t id = 0; id < num_variants; ++id) {
        string variant = ref_seq;

        // 1. Apply probabilistic substitutions per position
        for (size_t p = 0; p < variant.size(); ++p) {
            int roll = prob_dist(rng);
            if ((is_variable[p] && roll <= static_cast<int>(var_prob)) ||
                (!is_variable[p] && roll <= static_cast<int>(fw_prob))) {
                variant[p] = random_base(variant[p], rng);
            }
        }

        // 2. Fluctuate length at both ends (trim or extend with random bases)
        int delta_start = len_flunc_dist(rng);
        int delta_end   = len_flunc_dist(rng);
        size_t trim_start = delta_start > 0 ? min<size_t>(delta_start, variant.size()) : 0;
        size_t trim_end   = delta_end   > 0 ? min<size_t>(delta_end,   variant.size() - trim_start) : 0;

        string core = variant.substr(trim_start, variant.size() - trim_start - trim_end);

        if (delta_start < 0)
            core = random_bases(static_cast<size_t>(-delta_start), rng) + core;
        if (delta_end < 0)
            core += random_bases(static_cast<size_t>(-delta_end), rng);

        // 3. (Optional) reverse‑complement
        string out_seq = revcomp_flag ? reverse_complement(core) : core;

        // 4. Prepend barcode flanking sequence
        string current_barcode;
        if (barcode_flag) {
            current_barcode = all_barcode_ids[id % all_barcode_ids.size()];
        }
        if (barcode_flag) {
            const auto &bp = BARCODE_SEQS.at(current_barcode);
            string tag_seq;
            if (revcomp_flag) {
                tag_seq = "GGTGCTG" + bp.second + "TTAACCTTAGCAAT";
            } else {
                tag_seq = "AAGGTTAA" + bp.first + "CAGCACCT";
            }
            out_seq = tag_seq + out_seq;
        }

        // Build header tag
        string tag = revcomp_flag ? "_RC" : "_FR";
        if (barcode_flag)
            tag += "_" + current_barcode;

        // 5. Write to output FASTA (IDs start from 1)
        out << '>' << ref_name << tag << '_' << (id + 1) << '\n';
        out << out_seq << '\n';
    }

    return 0;
}
