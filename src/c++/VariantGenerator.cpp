// VariantGenerator: produces synthetic FASTA amplicons used to sanity-check
// and stress-test the downstream basecalling / alignment / DMS-filtering
// pipeline. A reference sequence is perturbed at two independent per-base
// rates (inside vs. outside user-defined mutation regions) and optionally
// extended or trimmed at both ends by a random length delta. With --revcomp
// a matching reverse-complement variant is emitted alongside each forward
// variant so the alignment stage can be tested against either orientation.
//
// Usage: VariantGenerator <var_prob> <fw_prob> <max_len_flux> <n>
//                         <region_file> <ref.fasta> <out.fasta>
//                         [--revcomp|-r]
// Probabilities are integer percentages in [0, 100]. region_file is a CSV
// list of `start,end` intervals (end exclusive) flagging the "variable" mask.
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>
#include <random>
#include <algorithm>
#include <cctype>
#include <stdexcept>

using namespace std;

// Half-open reference interval: [start, end).
struct Region {
    size_t start;   // inclusive
    size_t end;     // exclusive (must satisfy start < end)
};

// ──────────────────────────────────────────────────────────────────────────────
// Utility helpers
// ──────────────────────────────────────────────────────────────────────────────
// Draw a uniform random DNA base different from `exclude` (case-insensitive).
// Rejection sampling; expected number of iterations is ~4/3 per call.
static char random_base(char exclude, mt19937 &rng) {
    static const char bases[] = {'A', 'C', 'G', 'T'};
    uniform_int_distribution<int> dist(0, 3);
    char b;
    do {
        b = bases[dist(rng)];
    } while (toupper(b) == toupper(exclude));
    return b;
}

// Emit an `n`-long uniform random DNA string (used to pad reads on either end
// when `delta_start`/`delta_end` are negative, i.e. the variant is longer
// than the reference).
static string random_bases(size_t n, mt19937 &rng) {
    static const char bases[] = {'A', 'C', 'G', 'T'};
    uniform_int_distribution<int> dist(0, 3);
    string s;
    s.reserve(n);
    for (size_t i = 0; i < n; ++i)
        s.push_back(bases[dist(rng)]);
    return s;
}

// Watson-Crick complement of a single base; anything non-ACGT collapses to 'N'.
static inline char comp(char b) {
    switch (toupper(b)) {
        case 'A': return 'T';
        case 'C': return 'G';
        case 'G': return 'C';
        case 'T': return 'A';
        default:  return 'N';
    }
}

// Reverse-complement a DNA string by walking the input backwards and
// complementing each base. Keeps output length == input length.
static string reverse_complement(const string &seq) {
    string rc;
    rc.reserve(seq.size());
    for (auto it = seq.rbegin(); it != seq.rend(); ++it)
        rc.push_back(comp(*it));
    return rc;
}

// ──────────────────────────────────────────────────────────────────────────────
// FASTA & Region helpers
// ──────────────────────────────────────────────────────────────────────────────
// Read a single-record FASTA and return the concatenated sequence. Writes the
// contig name (the portion of the header after '>') into `sequence_name`.
// Throws on an unreadable file or an empty/nameless record.
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

// Parse a CSV-per-line "start,end" region file and flatten it into a vector
// of individual 0-based reference positions marked as variable. `end` is
// exclusive so [40,75) yields 35 positions. Throws on any interval that
// violates start<end or overruns the reference.
static vector<size_t> read_regions(const string &path, size_t ref_len) {
    ifstream in(path);
    if (!in)
        throw runtime_error("Cannot open region file: " + path);

    vector<size_t> allowed;
    string line;
    while (getline(in, line)) {
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
    //   8  [--revcomp|-r] (optional) – when specified, an *additional* set of
    //      reverse‑complement variants is generated, effectively doubling the
    //      output size.

    bool revcomp_flag = false;

    if (argc < 8) {
        cerr << "Usage: " << argv[0]
             << " <variable_mutant_prob> <framework_mutant_prob> <max_length_flunctuation>"
                " <num_variants> <region_file> <reference_fasta> <output_fasta>"
                " [--revcomp|-r]\n";
        return 1;
    }

    // Positional arguments
    const size_t var_prob        = stoul(argv[1]);
    const size_t fw_prob         = stoul(argv[2]);
    const size_t max_len_flunc   = stoul(argv[3]);
    size_t num_variants    = stoul(argv[4]);
    const string region_file     = argv[5];
    const string ref_fasta       = argv[6];
    const string out_fasta       = argv[7];

    // Optional flags
    for (int i = 8; i < argc; ++i) {
        string flag = argv[i];
        if (flag == "--revcomp" || flag == "-r") {
            revcomp_flag = true;
        } else {
            cerr << "Unknown flag: " << flag << '\n';
            return 1;
        }
    }

    // Validate probabilities
    if (var_prob > 100 || fw_prob > 100) {
        cerr << "Mutation probabilities must be in the range 0‑100\n";
        return 1;
    }

    // Reference sequence
    string ref_name;
    string ref_seq = read_fasta(ref_fasta, ref_name);
    const size_t ref_len = ref_seq.size();

    // Flatten region intervals into a per-base boolean mask for O(1) lookup
    // inside the hot mutation loop.
    vector<size_t> allowed = read_regions(region_file, ref_len);
    if (allowed.empty()) {
        cerr << "No regions specified\n";
        return 1;
    }
    vector<char> is_variable(ref_len, 0);
    for (size_t p : allowed) is_variable[p] = 1;

    // RNG helpers. `prob_dist` draws in [1,100] so that a probability threshold
    // of p means roll<=p, exactly matching the "p percent of the time" reading.
    // `len_flunc_dist` is signed so positive deltas trim and negative deltas
    // extend the ends.
    mt19937 rng(random_device{}());
    uniform_int_distribution<int> prob_dist(1, 100);  // inclusive
    uniform_int_distribution<int> len_flunc_dist(-static_cast<int>(max_len_flunc), static_cast<int>(max_len_flunc));

    // Build one variant by (1) point-mutating each base using the zone-specific
    // rate, (2) independently trimming or padding each end by a signed delta,
    // (3) optionally reverse-complementing the whole thing.
    auto generate_variant = [&](bool make_revcomp) -> string {
        string variant = ref_seq;

        // 1. Point mutations: inside-zone positions use `var_prob`, outside use
        //    the lower `fw_prob`. random_base guarantees a true substitution
        //    (new base != old base).
        for (size_t p = 0; p < variant.size(); ++p) {
            int roll = prob_dist(rng);
            if ((is_variable[p] && roll <= static_cast<int>(var_prob)) ||
                (!is_variable[p] && roll <= static_cast<int>(fw_prob))) {
                variant[p] = random_base(variant[p], rng);
            }
        }

        // 2. Length fluctuations at both ends. Positive delta -> trim that
        //    many bases from the corresponding end; negative delta -> pad with
        //    random bases. trim_end is clamped against the remaining length
        //    after start-trim to avoid underflow on the substr computation.
        int delta_start = len_flunc_dist(rng);
        int delta_end   = len_flunc_dist(rng);
        size_t trim_start = delta_start > 0 ? min<size_t>(delta_start, variant.size()) : 0;
        size_t trim_end   = delta_end   > 0 ? min<size_t>(delta_end,   variant.size() - trim_start) : 0;

        string core = variant.substr(trim_start, variant.size() - trim_start - trim_end);

        if (delta_start < 0)
            core = random_bases(static_cast<size_t>(-delta_start), rng) + core;
        if (delta_end < 0)
            core += random_bases(static_cast<size_t>(-delta_end), rng);

        // 3. Reverse-complement if requested for this specific variant.
        return make_revcomp ? reverse_complement(core) : core;
    };

    ofstream out(out_fasta);
    if (!out) {
        cerr << "Cannot open output file: " << out_fasta << "\n";
        return 1;
    }

    // ------------------------------------------------------------------
    // Variant generation loop
    // ------------------------------------------------------------------
    // With --revcomp on, pair every forward variant with an independent RC
    // variant (not the RC of the same draw); this doubles effective coverage
    // and exercises both strand orientations through the alignment stage.
    if(revcomp_flag){
        num_variants *= 2;  // Double the number of variants if reverse-complement is requested
        for(size_t id = 0; id < num_variants; id+=2){
            // Forward variant
            string fwd_seq = generate_variant(false);
            out << '>' << ref_name << "_FR_" << (id + 1) << '\n';
            out << fwd_seq << '\n';

            // Reverse-complement variant (independent draw, not RC of fwd_seq)
            string rc_seq = generate_variant(true);
            out << '>' << ref_name << "_RC_" << (id + 2) << '\n';
            out << rc_seq << '\n';
        }
    } else {
        for(size_t id = 0; id < num_variants; ++id) {
            // Forward variant
            string fwd_seq = generate_variant(false);
            out << '>' << ref_name << "_FR_" << (id + 1) << '\n';
            out << fwd_seq << '\n';
        }
    }

    return 0;
}
