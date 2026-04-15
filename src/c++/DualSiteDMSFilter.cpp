// DualSiteDMSFilter: flagship thesis tool.
//
// Produces the `DualSiteDMSFilter` executable. Reads a minimap2/htslib BAM of
// nanopore basecalls aligned to a single-contig amplicon reference, walks each
// read's CIGAR to reconstruct a reference-length "corrected" sequence, and
// emits two outputs:
//
//   1. A FASTQ (<out.fastq>) of reference-length polished reads. Bases inside
//      the user-supplied DMS zones are taken from the read (iff they pass
//      --min-bq); bases outside those zones are forced to the reference
//      ("WT") and stamped with the --wt-q phred quality. This is the input
//      consumed by the Benchmarker family of tools for accuracy evaluation.
//   2. Optionally, a TSV (--out-counts) of variant_key -> count pairs. The
//      variant_key is a sorted, ';'-joined list of "<pos><ref>>::<alt>" edits
//      within DMS zones ("WT" when no passing edits exist); this is the input
//      consumed by VariantConcordance.
//
// Core filters (in order): mapping flags (drop unmapped/secondary/supplementary
// and --min-mapq), CIGAR indel budget, zone coverage (full coverage of first+
// last DMS zone, or --require-all-zones for every zone), at most two mutated
// DMS codons per read ("dual-site" rule), and optional NNK codon validation
// (3rd base must be G or T; see --mut-codon-library NNK).
//
// Also accumulates per-region mismatch rates (DMS vs. framework / WT noise),
// per-base and per-codon, for stderr reporting. Framework noise is measured
// from the raw aligned bases outside DMS zones (not from the corrected_seq,
// which is forced WT there).
#include <htslib/sam.h>
#include <htslib/faidx.h>
#include <algorithm>
#include <cctype>
#include <cstdint>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>
#include <iomanip>

// 1-based inclusive reference interval describing one DMS mutation window.
// The thesis amplicon uses four zones (40-75, 211-243, 439-459, 577-627).
struct Zone {
    int32_t start; // 1-based inclusive
    int32_t end;   // 1-based inclusive
};

// Codon library design constraint applied to every mutated codon inside DMS
// zones. NNN accepts any codon; NNK (the thesis phage-display library design)
// additionally requires the 3rd base to be G or T.
enum class MutantCodonLibrary {
    NNN,
    NNK
};

// Convert all nucleotides to upper case. In-place via the value parameter; the
// unsigned-char cast avoids UB on signed chars with the high bit set.
static inline std::string to_upper(std::string s) {
    for (char &c : s) c = (char)std::toupper((unsigned char)c);
    return s;
}

// Parse the continuous zones literal (e.g., 11-40,50-80,90-130) into a vector of zones.
// Zones are sorted by start for deterministic first_zone/last_zone semantics.
// Returns false on any malformed token (missing dash, non-positive bounds,
// reversed range, or empty input) without touching `zones` beyond clearing it.
static bool parse_zones(const std::string& s, std::vector<Zone>& zones) {
    zones.clear();
    std::stringstream ss(s);
    std::string token;
    while (std::getline(ss, token, ',')) {
        auto dash = token.find('-');
        if (dash == std::string::npos) return false;
        std::string a = token.substr(0, dash);
        std::string b = token.substr(dash + 1);
        if (a.empty() || b.empty()) return false;
        int32_t st = std::stoi(a);
        int32_t en = std::stoi(b);
        if (st <= 0 || en <= 0 || en < st) return false;
        zones.push_back({st, en});
    }
    if (zones.empty()) return false;
    std::sort(zones.begin(), zones.end(), [](const Zone& x, const Zone& y){
        return x.start < y.start;
    });
    return true;
}

// Decode the 4-bit packed base at 0-based read offset `read_i0` from a BAM
// record. The `nt` lookup is htslib's standard IUPAC table indexed by the
// nibble value returned by bam_seqi.
static inline char bam_base(const bam1_t* b, int32_t read_i0) {
    static const char* nt = "=ACMGRSVTWYHKDBN";
    uint8_t* seq = bam_get_seq(b);
    char c = nt[bam_seqi(seq, read_i0)];
    return (char)std::toupper((unsigned char)c);
}

// Raw (non phred33-shifted) base quality at read offset `read_i0`. 255 is the
// htslib sentinel for "no quality recorded"; we coerce to 0 so downstream
// thresholds treat it as the lowest possible quality.
static inline uint8_t bam_qual_raw(const bam1_t* b, int32_t read_i0) {
    uint8_t q = bam_get_qual(b)[read_i0];
    if (q == 255) q = 0; // treat missing as 0
    return q;
}

// Join a sorted list of "<pos><ref>><alt>" edit tokens into the variant key
// format consumed by VariantConcordance. An empty edit list (no DMS-zone
// deviations after filtering) is encoded as the literal string "WT".
static std::string join_edits(const std::vector<std::string>& edits) {
    if (edits.empty()) return "WT";
    std::string out;
    out.reserve(64);
    for (size_t i = 0; i < edits.size(); i++) {
        if (i) out.push_back(';');
        out += edits[i];
    }
    return out;
}

// Parse the --mut-codon-library string into the enum. Returns false on
// anything other than "NNN" or "NNK" (case-insensitive).
static bool parse_mutant_codon_library(const std::string& s, MutantCodonLibrary& out) {
    std::string u = to_upper(s);
    if (u == "NNN") {
        out = MutantCodonLibrary::NNN;
        return true;
    }
    if (u == "NNK") {
        out = MutantCodonLibrary::NNK;
        return true;
    }
    return false;
}

// Pretty-printer used only in the stderr summary; keeps the log lines
// self-describing without needing to remember the enum value.
static const char* mutant_codon_library_name(MutantCodonLibrary lib) {
    switch (lib) {
        case MutantCodonLibrary::NNN: return "NNN";
        case MutantCodonLibrary::NNK: return "NNK";
    }
    return "unknown";
}

// Print the full usage/help string including all optional flags. Using a raw
// string literal keeps quoting tractable for a multi-line block.
static void usage() {
    std::cerr <<
R"(Usage:
  DualSiteDMSFilter <in.bam> <ref.fasta> <out.fastq> <zones>
    --codon-frame <0|1|2>            (default: 0; 0 means ref pos 1 starts a codon)
    --max-indel-events <N>           (default: 0)
    --max-indel-bases  <N>           (default: 0)
    --min-mapq <Q>                   (default: 0)
    --min-bq <Q>                     (default: 0; base quality threshold inside DMS zones)
    --min-fw-bq <Q>                  (default: same as --min-bq; base quality threshold for framework noise)
    --wt-q <Q>                       (default: 40; phred for WT-filled positions in corrected FASTQ)
    --mut-codon-library <NNN|NNK>    (default: NNN; apply library design rule to mutated codons)
    --require-all-zones              (default: off; if set, require full coverage of ALL zones, not just first+last)
    --out-counts <tsv>               (optional: writes formatted variant key)

Notes:
  - corrected FASTQ sequences are ALWAYS reference-length.
  - outside DMS zones = WT (reference), even if read disagrees.
  - qualities for WT-filled bases are set to --wt-q.
  - dual-site DMS filter: reads may carry at most 2 mutated codons across the DMS zones.
  - if --mut-codon-library NNK is used, every mutated codon must be fully observed at
    passing base quality and have a 3rd base of G or T.
  - framework (WT) noise is measured from RAW aligned read bases vs reference OUTSIDE zones,
    not from corrected_seq (which is forced WT outside zones).
  - "per-base mismatch rates" are computed as mismatches / bases_compared for each region.
  - "per-codon mismatch rates" are computed as codon_mismatches / codons_compared for each region,
    where a codon is "compared" if at least one base in that codon passes filters in that region.
)";
}

// Return true iff a single interval in `merged` fully covers zone `z`. `merged`
// must already be sorted and gap-merged. Terminates on the first interval that
// could possibly contain `z.start`: if that interval starts after `z.start`
// there is a gap, otherwise the decision reduces to end >= z.end.
static bool zone_fully_covered(const Zone& z, const std::vector<std::pair<int32_t,int32_t>>& merged) {
    for (const auto& iv : merged) {
        if (iv.second < z.start) continue;
        if (iv.first > z.start) return false;
        return iv.second >= z.end;
    }
    return false;
}

// Encode a raw phred quality value into the ASCII phred+33 FASTQ representation.
// Clamps at 93 so the output stays printable ('!'..'~').
static inline char phred33(uint8_t q) {
    if (q > 93) q = 93;
    return (char)(q + 33);
}

// Safe codon index: avoids negative integer division weirdness.
// Returns -1 if (pos-1) is before the first codon start based on codon_frame.
// Guarding against the negative-dividend case keeps the function safe even
// when the reference position falls in the short pre-frame prefix (e.g.
// codon_frame=2 means positions 1-2 have no codon).
static inline int32_t codon_index_0based(int32_t pos_1based, int codon_frame, int32_t n_codons) {
    int32_t x = pos_1based - 1; // 0-based base index
    if (x < codon_frame) return -1;
    int32_t idx = (x - codon_frame) / 3;
    if (idx < 0 || idx >= n_codons) return -1;
    return idx;
}

// Inverse of codon_index_0based: returns the 1-based reference position of the
// first base in codon `codon_idx_0based` under the given frame.
static inline int32_t codon_start_1based(int32_t codon_idx_0based, int codon_frame) {
    return codon_frame + 1 + (codon_idx_0based * 3);
}

// True iff `c` is one of the four canonical DNA bases (uppercase only; caller
// is expected to have already folded case).
static inline bool is_acgt(char c) {
    return c == 'A' || c == 'C' || c == 'G' || c == 'T';
}

// Validate that a mutated codon conforms to the configured library design.
// For NNN every codon is accepted. For NNK every one of the three codon bases
// must (a) live inside a DMS zone, (b) have been observed in the read at
// passing base quality (zone_hq_seen), (c) be a canonical ACGT base, and
// (d) the 3rd base must be G or T (the defining NNK constraint). Any violated
// condition returns false so the caller can reject the whole read.
static bool mutated_codon_matches_library(
    const std::string& corrected_seq,
    const std::vector<uint8_t>& zone_mask,
    const std::vector<uint8_t>& zone_hq_seen,
    int32_t codon_start,
    int32_t ref_len,
    MutantCodonLibrary library
) {
    if (library == MutantCodonLibrary::NNN) return true;
    if (codon_start < 1 || codon_start + 2 > ref_len) return false;

    // All three positions must be inside a DMS zone AND observed at HQ;
    // otherwise we cannot trust the 3rd-base constraint at this codon.
    for (int32_t p = codon_start; p <= codon_start + 2; ++p) {
        if (!zone_mask[p - 1] || !zone_hq_seen[p - 1]) return false;
    }

    char b1 = corrected_seq[codon_start - 1];
    char b2 = corrected_seq[codon_start];
    char b3 = corrected_seq[codon_start + 1];
    if (!is_acgt(b1) || !is_acgt(b2) || !is_acgt(b3)) return false;

    // NNK constraint: 3rd base must be G or T.
    return b3 == 'G' || b3 == 'T';
}

int main(int argc, char** argv) {
    if (argc < 5) { usage(); return 3; }

    std::string in_bam = argv[1];
    std::string ref_fa = argv[2];
    std::string out_fq = argv[3];

    std::vector<Zone> zones;
    if(!parse_zones(argv[4], zones)){
        std::cerr << "Failed to parse zones: " << argv[4] << "\n";
        return 2;
    }
    if (zones.empty()) {
        std::cerr << "Error: zones cannot be empty.\n";
        usage();
        return 2;
    }

    // Defaults for every optional flag. -1 on min_fw_bq is a sentinel meaning
    // "fall back to min_bq if the user doesn't override it"; wt_q=40 gives
    // WT-filled bases a plausible Illumina-quality-like phred so downstream
    // tools do not mistake them for sequencing noise.
    int codon_frame = 0;
    int max_indel_events = 0;
    int max_indel_bases  = 0;
    int min_mapq = 0;
    int min_bq   = 0;
    int min_fw_bq = -1; // default: same as min_bq
    int wt_q     = 40;
    MutantCodonLibrary mutant_codon_library = MutantCodonLibrary::NNN;
    bool require_all_zones = false;
    std::string out_counts;

    for (int i = 5; i < argc; i++) {
        std::string a = argv[i];
        auto need = [&](const char* opt) -> const char* {
            if (i + 1 >= argc) {
                std::cerr << "Missing value for " << opt << "\n";
                std::exit(2);
            }
            return argv[++i];
        };
        if (a == "--codon-frame") codon_frame = std::stoi(need("--codon-frame"));
        else if (a == "--max-indel-events") max_indel_events = std::stoi(need("--max-indel-events"));
        else if (a == "--max-indel-bases")  max_indel_bases  = std::stoi(need("--max-indel-bases"));
        else if (a == "--min-mapq")         min_mapq         = std::stoi(need("--min-mapq"));
        else if (a == "--min-bq")           min_bq           = std::stoi(need("--min-bq"));
        else if (a == "--min-fw-bq")        min_fw_bq        = std::stoi(need("--min-fw-bq"));
        else if (a == "--wt-q")             wt_q             = std::stoi(need("--wt-q"));
        else if (a == "--mut-codon-library") {
            if (!parse_mutant_codon_library(need("--mut-codon-library"), mutant_codon_library)) {
                std::cerr << "Unsupported --mut-codon-library value. Use NNN or NNK.\n";
                return 2;
            }
        }
        else if (a == "--require-all-zones") require_all_zones = true;
        else if (a == "--out-counts") out_counts = need("--out-counts");
        else {
            std::cerr << "Unknown option: " << a << "\n";
            usage();
            return 2;
        }
    }

    // Frame must live in {0,1,2}; silently clamp rather than error so that a
    // downstream wrapper script passing an arithmetic expression can not crash
    // the pipeline on a minor miscalculation.
    if (codon_frame < 0) codon_frame = 0;
    if (codon_frame > 2) codon_frame = 2;

    if (min_fw_bq < 0) min_fw_bq = min_bq;

    if (wt_q < 0) wt_q = 0;
    if (wt_q > 93) wt_q = 93;

    // Load reference FASTA: require exactly ONE contig, and use it.
    // fai_load auto-generates or refreshes the .fai sidecar if missing.
    faidx_t* fai = fai_load(ref_fa.c_str());
    if (!fai) {
        std::cerr << "Failed to load reference FASTA index. Run: samtools faidx " << ref_fa << "\n";
        return 2;
    }

    // Single-contig constraint. The tool assumes one amplicon reference and
    // uses zone coordinates relative to it; multi-contig input would make the
    // zone-bounds check ambiguous.
    int nseq = faidx_nseq(fai);
    if (nseq != 1) {
        std::cerr << "Reference FASTA must contain exactly 1 contig, but found " << nseq << ".\n";
        std::cerr << "Fix: provide a single-contig FASTA (your amplicon reference only).\n";
        fai_destroy(fai);
        return 2;
    }

    const char* contig = faidx_iseq(fai, 0);
    if (!contig || !*contig) {
        std::cerr << "Failed to read contig name from FASTA index.\n";
        fai_destroy(fai);
        return 2;
    }

    // htslib returns int64_t lengths; the DMS logic below uses int32_t for
    // positions so reject any reference longer than INT32_MAX up front rather
    // than risk silent overflow in codon math.
    int64_t ref_len64 = faidx_seq_len(fai, contig);
    if (ref_len64 <= 0) {
        std::cerr << "Failed to get reference length for contig: " << contig << "\n";
        fai_destroy(fai);
        return 2;
    }
    if (ref_len64 > INT32_MAX) {
        std::cerr << "Reference contig is too long for this tool (len=" << ref_len64 << ").\n";
        fai_destroy(fai);
        return 2;
    }
    int32_t ref_len = (int32_t)ref_len64;

    int fetched_len = 0;
    char* ref_c = faidx_fetch_seq(fai, contig, 0, ref_len - 1, &fetched_len);
    if (!ref_c || fetched_len != ref_len) {
        std::cerr << "Failed to fetch full reference sequence for " << contig
                  << " expected len=" << ref_len << " got=" << fetched_len << "\n";
        if (ref_c) free(ref_c);
        fai_destroy(fai);
        return 2;
    }

    std::string ref_full = to_upper(std::string(ref_c, ref_c + fetched_len));
    free(ref_c);

    // 0-based boolean mask: zone_mask[p] == 1 iff reference position p+1
    // lives inside any DMS zone. Flattening here keeps the hot CIGAR loop's
    // zone test to an O(1) array lookup per base.
    std::vector<uint8_t> zone_mask(ref_len, 0);
    for (const auto& z : zones) {
        if (z.start < 1 || z.end > ref_len) {
            std::cerr << "Zone out of bounds for reference length " << ref_len
                      << ": " << z.start << "-" << z.end << "\n";
            fai_destroy(fai);
            return 2;
        }
        for (int32_t p = z.start; p <= z.end; p++) zone_mask[p - 1] = 1;
    }

    // zones is already sorted in parse_zones, so first/last are by position.
    // The default coverage rule only requires both endpoints, not every zone.
    Zone first_zone = zones.front();
    Zone last_zone  = zones.back();

    // Precompute zone positions (1-based) once. Enumerating zone_pos avoids
    // re-scanning zone_mask for every read during variant-key construction.
    std::vector<int32_t> zone_pos;
    zone_pos.reserve(ref_len);
    for (int32_t p = 1; p <= ref_len; ++p) {
        if (zone_mask[p - 1]) zone_pos.push_back(p);
    }

    // Open BAM
    samFile* in = sam_open(in_bam.c_str(), "r");
    if (!in) {
        std::cerr << "Failed to open BAM: " << in_bam << "\n";
        fai_destroy(fai);
        return 2;
    }
    bam_hdr_t* hdr = sam_hdr_read(in);
    if (!hdr) {
        std::cerr << "Failed to read BAM header.\n";
        sam_close(in);
        fai_destroy(fai);
        return 2;
    }

    std::ofstream out_fq_stream(out_fq);
    if (!out_fq_stream) {
        std::cerr << "Failed to open corrected FASTQ: " << out_fq << "\n";
        bam_hdr_destroy(hdr);
        sam_close(in);
        fai_destroy(fai);
        return 2;
    }

    int contig_tid = bam_name2id(hdr, contig);
    if (contig_tid < 0) {
        std::cerr << "Contig not found in BAM header: " << contig << "\n";
        out_fq_stream.close();
        bam_hdr_destroy(hdr);
        sam_close(in);
        fai_destroy(fai);
        return 2;
    }

    bam1_t* b = bam_init1();
    if (!b) {
        std::cerr << "Failed to allocate bam record.\n";
        out_fq_stream.close();
        bam_hdr_destroy(hdr);
        sam_close(in);
        fai_destroy(fai);
        return 2;
    }

    // If --out-counts is set we accumulate per-variant_key counts; otherwise
    // we only track which keys appeared, purely to report `unique_variants` in
    // the stderr summary. Reserve 200k buckets to front-load allocation on a
    // typical thesis-scale dataset (~100k-200k distinct DMS variants).
    const bool want_table = !out_counts.empty();
    std::unordered_map<std::string, uint64_t> counts;
    std::unordered_set<std::string> unique_keys;
    if (want_table) counts.reserve(200000);
    else unique_keys.reserve(200000);

    // Stats
    uint64_t n_total = 0;
    uint64_t n_skip_tid = 0;
    uint64_t n_fail_mapq = 0;
    uint64_t n_fail_indel = 0;
    uint64_t n_fail_coverage = 0;
    uint64_t n_fail_dual_site = 0;
    uint64_t n_fail_mutant_codon_library = 0;
    uint64_t n_pass = 0;
    uint64_t n_pass_0_codon_mut = 0;
    uint64_t n_pass_1_codon_mut = 0;
    uint64_t n_pass_2_codon_mut = 0;

    // Zone mismatch sums (signal, inside mutation regions) - per read sums
    uint64_t sum_nt_mismatches = 0;
    uint64_t sum_codon_mismatches = 0;

    // Framework mismatch sums (noise, outside mutation regions) - per read sums
    uint64_t sum_fw_nt_mismatches = 0;
    uint64_t sum_fw_codon_mismatches = 0;

    // NEW: denominators for mismatch rates (PASS reads only)
    uint64_t sum_zone_bases_compared = 0;
    uint64_t sum_fw_bases_compared   = 0;
    uint64_t sum_zone_codons_compared = 0;
    uint64_t sum_fw_codons_compared   = 0;

    // Number of complete codons fitting in the reference (ceil division).
    const int32_t n_codons = (ref_len + 2) / 3;

    // "Epoch-stamped" visited arrays. Instead of clearing a size-n_codons
    // vector per read, we store a per-codon counter and compare against the
    // current epoch; a codon is "seen this read" iff slot == epoch. Bumping
    // epoch at the start of each read is O(1); we only refill the array when
    // epoch would wrap to 0. One pair of arrays per category: DMS-zone
    // mismatches, framework mismatches, DMS-zone codons-compared denominator,
    // framework codons-compared denominator.
    std::vector<uint32_t> codon_seen_mis((size_t)n_codons, 0);
    uint32_t epoch_mis = 1;

    std::vector<uint32_t> codon_seen_fw_mis((size_t)n_codons, 0);
    uint32_t epoch_fw_mis = 1;

    // NEW: Codon markers for "codons compared" (zone and framework)
    std::vector<uint32_t> codon_seen_cmp((size_t)n_codons, 0);
    uint32_t epoch_cmp = 1;

    std::vector<uint32_t> codon_seen_fw_cmp((size_t)n_codons, 0);
    uint32_t epoch_fw_cmp = 1;

    while (sam_read1(in, hdr, b) >= 0) {
        n_total++;

        // Primary-alignment gate: drop unmapped / secondary / supplementary
        // records. Equivalent to samtools view -F 0x904. These records either
        // carry no usable sequence or duplicate bases already counted against
        // the primary alignment.
        if (b->core.flag & (BAM_FUNMAP | BAM_FSECONDARY | BAM_FSUPPLEMENTARY)) continue;
        if (b->core.tid != contig_tid) { n_skip_tid++; continue; }
        if ((int)b->core.qual < min_mapq) { n_fail_mapq++; continue; }

        uint32_t* cigar = bam_get_cigar(b);

        // INDEL assessment. Walk the CIGAR once to count insertion+deletion
        // events and bases, and reject the read if either exceeds its budget.
        // The inner break exits early once the budgets are exceeded; the
        // second check is then just a formality that still returns the
        // early-exit state correctly.
        int indel_events = 0;
        int indel_bases = 0;
        for (uint32_t i = 0; i < b->core.n_cigar; i++) {
            int op = bam_cigar_op(cigar[i]);
            int len = bam_cigar_oplen(cigar[i]);
            if (op == BAM_CINS || op == BAM_CDEL) {
                indel_events++;
                indel_bases += len;
                if (indel_events > max_indel_events || indel_bases > max_indel_bases) break;
            }
        }
        if (indel_events > max_indel_events || indel_bases > max_indel_bases) {
            n_fail_indel++;
            continue;
        }

        // Start CIGAR walk. refpos tracks the 1-based reference coordinate,
        // readpos the 0-based offset into the packed read sequence. `cov`
        // collects M/=/X intervals so we can later check zone coverage.
        int32_t refpos = b->core.pos + 1; // 1-based
        int32_t readpos = 0;

        std::vector<std::pair<int32_t,int32_t>> cov;
        cov.reserve(b->core.n_cigar);

        // Pre-fill the corrected FASTQ with the reference and a uniform
        // WT quality string. Only DMS-zone positions with HQ read bases will
        // overwrite these defaults below; outside zones the reference always
        // wins, which is exactly the DMS-aware polishing contract.
        std::string corrected_seq = ref_full; // WT baseline
        std::string corrected_qual;
        corrected_qual.assign(ref_len, phred33((uint8_t)wt_q));
        std::vector<uint8_t> zone_hq_seen((size_t)ref_len, 0);

        // Per-read framework mismatch counters (raw read base vs ref outside zones)
        uint32_t read_fw_nt_mis = 0;
        uint32_t read_fw_codon_mis = 0;

        // NEW: per-read denominators
        uint32_t read_zone_bases_cmp = 0;
        uint32_t read_fw_bases_cmp   = 0;
        uint32_t read_zone_codons_cmp = 0;
        uint32_t read_fw_codons_cmp   = 0;

        // Assign a fresh epoch id per read for each codon-tracking category.
        // On wrap-around (epoch hits 0 after increment) the slots must be
        // zeroed so the next epoch's comparisons start from a clean state.
        uint32_t this_epoch_fw_mis = epoch_fw_mis++;
        if (epoch_fw_mis == 0) {
            std::fill(codon_seen_fw_mis.begin(), codon_seen_fw_mis.end(), 0);
            epoch_fw_mis = 1;
            this_epoch_fw_mis = epoch_fw_mis++;
        }

        // epochs for codons compared
        uint32_t this_epoch_cmp = epoch_cmp++;
        if (epoch_cmp == 0) {
            std::fill(codon_seen_cmp.begin(), codon_seen_cmp.end(), 0);
            epoch_cmp = 1;
            this_epoch_cmp = epoch_cmp++;
        }

        uint32_t this_epoch_fw_cmp = epoch_fw_cmp++;
        if (epoch_fw_cmp == 0) {
            std::fill(codon_seen_fw_cmp.begin(), codon_seen_fw_cmp.end(), 0);
            epoch_fw_cmp = 1;
            this_epoch_fw_cmp = epoch_fw_cmp++;
        }

        // Walk CIGAR. Each op advances refpos/readpos according to SAM spec:
        //   M/=/X consume both; record coverage and inspect every aligned base
        //   I/S    consume read only  (insertion / soft clip)
        //   D/N    consume ref only   (deletion / ref skip)
        //   H/P    consume neither    (hard clip / pad)
        for (uint32_t ci = 0; ci < b->core.n_cigar; ci++) {
            int op  = bam_cigar_op(cigar[ci]);
            int len = bam_cigar_oplen(cigar[ci]);

            if (op == BAM_CMATCH || op == BAM_CEQUAL || op == BAM_CDIFF) {
                cov.push_back({refpos, refpos + len - 1});

                for (int k = 0; k < len; k++) {
                    int32_t rp = refpos + k;   // 1-based ref
                    int32_t qp = readpos + k;  // 0-based read
                    if (rp < 1 || rp > ref_len) continue;
                    if (qp < 0 || qp >= (int32_t)b->core.l_qseq) continue;

                    char rb = bam_base(b, qp);
                    uint8_t rq = bam_qual_raw(b, qp);

                    char refb = ref_full[rp - 1];
                    if (refb == 'N') {
                        // not comparable in either region
                        continue;
                    }

                    // Framework region: outside DMS zones. We measure noise
                    // here against the RAW read base (rb) vs. reference, not
                    // against corrected_seq (which is forced WT outside zones,
                    // so would trivially report zero framework mismatches).
                    if (!zone_mask[rp - 1]) {
                        if (rb != 'N' && (int)rq >= min_fw_bq) {
                            // denominator: bases compared
                            read_fw_bases_cmp++;

                            // denominator: codons compared (once per read per codon).
                            // The epoch check ensures a single codon contributes
                            // exactly one unit to the denominator even if all
                            // three of its bases pass the filter.
                            int32_t codon_idx = codon_index_0based(rp, codon_frame, n_codons);
                            if (codon_idx >= 0) {
                                uint32_t &slot_cmp = codon_seen_fw_cmp[(size_t)codon_idx];
                                if (slot_cmp != this_epoch_fw_cmp) {
                                    slot_cmp = this_epoch_fw_cmp;
                                    read_fw_codons_cmp++;
                                }
                            }

                            // mismatch against the reference base (noise signal)
                            if (rb != refb) {
                                read_fw_nt_mis++;

                                // codon mismatch (once per codon per read) —
                                // same epoch pattern as the denominator.
                                if (codon_idx >= 0) {
                                    uint32_t &slot_mis = codon_seen_fw_mis[(size_t)codon_idx];
                                    if (slot_mis != this_epoch_fw_mis) {
                                        slot_mis = this_epoch_fw_mis;
                                        read_fw_codon_mis++;
                                    }
                                }
                            }
                        }
                    }

                    // Mutational region: inside DMS zones. Here we BOTH
                    // accumulate the per-read zone denominators AND transfer
                    // the HQ read base into corrected_seq, flagging the
                    // position in zone_hq_seen so NNK codon validation later
                    // knows this position was observed at passing quality.
                    if (zone_mask[rp - 1]) {
                        if (rb != 'N' && (int)rq >= min_bq) {
                            // denominator: bases compared
                            read_zone_bases_cmp++;

                            // denominator: codons compared (once per read per codon)
                            int32_t codon_idx = codon_index_0based(rp, codon_frame, n_codons);
                            if (codon_idx >= 0) {
                                uint32_t &slot_cmp = codon_seen_cmp[(size_t)codon_idx];
                                if (slot_cmp != this_epoch_cmp) {
                                    slot_cmp = this_epoch_cmp;
                                    read_zone_codons_cmp++;
                                }
                            }

                            // Mutational region fill into corrected seq/qual (your existing logic)
                            corrected_seq[rp - 1]  = rb;
                            corrected_qual[rp - 1] = phred33(rq);
                            zone_hq_seen[rp - 1] = 1;
                        }
                    }
                }

                refpos += len;
                readpos += len;
            }
            else if (op == BAM_CINS || op == BAM_CSOFT_CLIP) {
                readpos += len;  // consumes read only
            }
            else if (op == BAM_CDEL || op == BAM_CREF_SKIP) {
                refpos += len;   // consumes ref only
            }
            else if (op == BAM_CHARD_CLIP || op == BAM_CPAD) {
                // consumes neither; no advancement
            }
        }

        // Merge M/=/X intervals into a minimal set of disjoint reference
        // spans so zone_fully_covered can short-circuit. Intervals that
        // touch (gap == 1 base, i.e. a.second+1 == b.first) are also merged,
        // which correctly treats a single-base deletion within coverage as
        // "still covered on both sides".
        if (cov.empty()) { n_fail_coverage++; continue; }
        std::sort(cov.begin(), cov.end());
        std::vector<std::pair<int32_t,int32_t>> merged;
        merged.reserve(cov.size());
        for (const auto& iv : cov) {
            if (merged.empty() || iv.first > merged.back().second + 1) merged.push_back(iv);
            else merged.back().second = std::max(merged.back().second, iv.second);
        }

        // Coverage rule. Default: demand full coverage only of the first and
        // last zone (by position) — a read that spans the whole DMS region
        // necessarily spans the intermediate zones too for the thesis
        // amplicon. --require-all-zones is stricter and verifies each zone
        // independently.
        bool ok = true;
        if (require_all_zones) {
            for (const auto& z : zones) {
                if (!zone_fully_covered(z, merged)) { ok = false; break; }
            }
        } else {
            if (!zone_fully_covered(first_zone, merged)) ok = false;
            if (ok && !zone_fully_covered(last_zone, merged)) ok = false;
        }

        if (!ok) {
            n_fail_coverage++;
            continue;
        }

        // Variant key + mismatch stats (zone-aware, quality-filtered via corrected_seq).
        // Mismatches at this stage are measured on corrected_seq vs. reference,
        // so positions that failed the HQ filter left the reference baseline
        // unchanged and therefore cannot contribute a spurious mismatch.
        uint32_t this_epoch_mis = epoch_mis++;
        if (epoch_mis == 0) {
            std::fill(codon_seen_mis.begin(), codon_seen_mis.end(), 0);
            epoch_mis = 1;
            this_epoch_mis = epoch_mis++;
        }

        uint32_t read_nt_mis = 0;
        uint32_t read_codon_mis = 0;

        // Edit tokens are accumulated in 1-based reference order (we iterate
        // zone_pos ascending). Reserve 8 to cover the typical 0-2 mutations
        // plus headroom. mutated_codon_indices caches the distinct codon
        // indices so NNK validation below can skip re-scanning.
        std::vector<std::string> edits;
        edits.reserve(8);
        std::vector<int32_t> mutated_codon_indices;
        mutated_codon_indices.reserve(2);

        for (int32_t p : zone_pos) {
            char refb = ref_full[p - 1];
            char cb   = corrected_seq[p - 1];

            if (cb != 'N' && refb != 'N' && cb != refb) {
                read_nt_mis++;

                int32_t codon_idx = codon_index_0based(p, codon_frame, n_codons);
                if (codon_idx >= 0) {
                    uint32_t &slot = codon_seen_mis[(size_t)codon_idx];
                    if (slot != this_epoch_mis) {
                        slot = this_epoch_mis;
                        read_codon_mis++;
                        mutated_codon_indices.push_back(codon_idx);
                    }
                }

                std::string e;
                e.reserve(16);
                e += std::to_string(p);
                e.push_back(refb);
                e.push_back('>');
                e.push_back(cb);
                edits.push_back(std::move(e));
            }
        }

        // Dual-site DMS library: allow at most two mutated codons per passing read.
        // Above this budget the read is more likely a real multi-mutation or
        // nanopore noise than the single/dual site the library was designed
        // to produce; dropping it sharply reduces false-positive variants.
        if (read_codon_mis > 2) {
            n_fail_dual_site++;
            continue;
        }

        // Validate every mutated codon against the configured library design.
        // For NNK this enforces the "3rd base ∈ {G,T}" constraint on top of
        // full HQ observation of all three codon bases.
        bool mutant_codons_ok = true;
        for (int32_t codon_idx : mutated_codon_indices) {
            int32_t codon_start = codon_start_1based(codon_idx, codon_frame);
            if (!mutated_codon_matches_library(
                    corrected_seq,
                    zone_mask,
                    zone_hq_seen,
                    codon_start,
                    ref_len,
                    mutant_codon_library)) {
                mutant_codons_ok = false;
                break;
            }
        }
        if (!mutant_codons_ok) {
            n_fail_mutant_codon_library++;
            continue;
        }

        // PASS read: framework noise, denominators, and zone mutation stats are valid to accumulate.
        sum_fw_nt_mismatches += read_fw_nt_mis;
        sum_fw_codon_mismatches += read_fw_codon_mis;

        sum_zone_bases_compared += read_zone_bases_cmp;
        sum_fw_bases_compared   += read_fw_bases_cmp;
        sum_zone_codons_compared += read_zone_codons_cmp;
        sum_fw_codons_compared   += read_fw_codons_cmp;

        sum_nt_mismatches += read_nt_mis;
        sum_codon_mismatches += read_codon_mis;

        if (read_codon_mis == 0) n_pass_0_codon_mut++;
        else if (read_codon_mis == 1) n_pass_1_codon_mut++;
        else n_pass_2_codon_mut++;

        // Sort the edits lexicographically before joining so that reads which
        // agree on the set of mutations generate identical variant_keys
        // regardless of the order zone_pos visited their positions.
        std::sort(edits.begin(), edits.end());
        std::string key = join_edits(edits);

        if (want_table) counts[key]++;
        else unique_keys.insert(std::move(key));

        // Write corrected FASTQ. Header preserves the original qname plus a
        // minimal set of provenance tags (mapq, contig) useful when eyeballing
        // the output. Sequence and quality are both exactly ref_len bytes.
        out_fq_stream << "@" << bam_get_qname(b)
                      << " mapq=" << (int)b->core.qual
                      << " contig=" << contig << "\n";
        out_fq_stream << corrected_seq << "\n+\n" << corrected_qual << "\n";

        n_pass++;
    }

    // Write counts TSV. Sorted by count descending, with lexicographic
    // variant_key as a tie-breaker, so the output is deterministic and
    // human-readable (most-abundant variants at top).
    if (!out_counts.empty()){
        std::vector<std::pair<std::string, uint64_t>> items;
        items.reserve(counts.size());
        for (const auto& kv : counts) items.push_back(kv);

        std::sort(items.begin(), items.end(), [](const auto& x, const auto& y){
            if (x.second != y.second) return x.second > y.second;
            return x.first < y.first;
        });

        std::ofstream ofs(out_counts);
        if (!ofs) {
            std::cerr << "Failed to open output counts file: " << out_counts << "\n";
            bam_destroy1(b);
            bam_hdr_destroy(hdr);
            sam_close(in);
            fai_destroy(fai);
            return 2;
        }
        ofs << "variant_key\tcount\n";
        for (const auto& it : items) ofs << it.first << "\t" << it.second << "\n";
    }

    std::cerr << std::setprecision(10);

    std::cerr << "[DmsFilterCount]\n";
    std::cerr << "in_bam=" << in_bam << "\n";
    std::cerr << "contig=" << contig << "\n";
    std::cerr << "ref_len=" << ref_len << "\n";
    std::cerr << "out_fastq=" << out_fq << "\n";
    std::cerr << "zones=";
    for (size_t i = 0; i < zones.size(); i++) {
        if (i) std::cerr << ",";
        std::cerr << zones[i].start << "-" << zones[i].end;
    }
    std::cerr << "\n";
    std::cerr << "require_all_zones=" << (require_all_zones ? "true" : "false") << "\n";
    std::cerr << "codon_frame=" << codon_frame << "\n";
    std::cerr << "mut_codon_library=" << mutant_codon_library_name(mutant_codon_library) << "\n";
    std::cerr << "max_indel_events=" << max_indel_events << " max_indel_bases=" << max_indel_bases << "\n";
    std::cerr << "min_mapq=" << min_mapq << " min_bq=" << min_bq << " min_fw_bq=" << min_fw_bq << " wt_q=" << wt_q << "\n";
    std::cerr << "total_reads=" << n_total << "\n";
    std::cerr << "skip_other_contig=" << n_skip_tid << "\n";
    std::cerr << "fail_mapq=" << n_fail_mapq << "\n";
    std::cerr << "fail_indel=" << n_fail_indel << "\n";
    std::cerr << "fail_coverage=" << n_fail_coverage << "\n";
    std::cerr << "fail_dual_site=" << n_fail_dual_site << "\n";
    std::cerr << "fail_mut_codon_library=" << n_fail_mutant_codon_library << "\n";
    std::cerr << "pass=" << n_pass << "\n";
    std::cerr << "pass_reads_0_codon_mutations=" << n_pass_0_codon_mut << "\n";
    std::cerr << "pass_reads_1_codon_mutations=" << n_pass_1_codon_mut << "\n";
    std::cerr << "pass_reads_2_codon_mutations=" << n_pass_2_codon_mut << "\n";
    std::cerr << "pass_rate=" << (n_total ? (static_cast<double>(n_pass) / n_total) : 0.0) << "\n";
    std::cerr << "unique_variants=" << (want_table ? counts.size() : unique_keys.size()) << "\n";

    // Per-pass-read averages. "Net" values subtract framework (noise) from
    // DMS (signal) and floor at zero; interpret as the average number of
    // DMS-specific mismatches per read above the background framework rate.
    double avg_nt_mis = (n_pass ? (double)sum_nt_mismatches / (double)n_pass : 0.0);
    double avg_codon_mis = (n_pass ? (double)sum_codon_mismatches / (double)n_pass : 0.0);
    double avg_fw_nt_mis = (n_pass ? (double)sum_fw_nt_mismatches / (double)n_pass : 0.0);
    double avg_fw_codon_mis = (n_pass ? (double)sum_fw_codon_mismatches / (double)n_pass : 0.0);
    // per-read net differences
    double net_avg_nt_mis = avg_nt_mis - avg_fw_nt_mis;
    double net_avg_codon_mis = avg_codon_mis - avg_fw_codon_mis;
    if (net_avg_nt_mis < 0) net_avg_nt_mis = 0.0;
    if (net_avg_codon_mis < 0) net_avg_codon_mis = 0.0;
    
    std::cerr << "avg_dms_nucleotide_mismatches_per_read=" << avg_nt_mis << "\n";
    std::cerr << "avg_framework_nucleotide_mismatches_per_read=" << avg_fw_nt_mis << "\n";
    std::cerr << "avg_net_nucleotide_mismatches_per_read=" << net_avg_nt_mis << "\n";
    std::cerr << "avg_dms_codon_mismatches_per_read=" << avg_codon_mis << "\n";
    std::cerr << "avg_framework_codon_mismatches_per_read=" << avg_fw_codon_mis << "\n";
    std::cerr << "avg_net_codon_mismatches_per_read=" << net_avg_codon_mis << "\n";

    // denominators
    std::cerr << "dms_bases_compared=" << sum_zone_bases_compared << "\n";
    std::cerr << "framework_bases_compared=" << sum_fw_bases_compared << "\n";
    std::cerr << "dms_codons_compared=" << sum_zone_codons_compared << "\n";
    std::cerr << "framework_codons_compared=" << sum_fw_codons_compared << "\n";

    // Per-base mismatch rates: mismatches / bases_compared for each region.
    // This is the normalised metric reported in the thesis figures and is
    // independent of read length, unlike the per-read averages above.
    double zone_nt_rate = (sum_zone_bases_compared ? (double)sum_nt_mismatches / (double)sum_zone_bases_compared : 0.0);
    double fw_nt_rate   = (sum_fw_bases_compared   ? (double)sum_fw_nt_mismatches / (double)sum_fw_bases_compared   : 0.0);
    double net_nt_rate  = zone_nt_rate - fw_nt_rate;
    if (net_nt_rate < 0) net_nt_rate = 0.0;

    std::cerr << "dms_nt_mismatch_rate_per_base=" << zone_nt_rate << "\n";
    std::cerr << "framework_nt_mismatch_rate_per_base=" << fw_nt_rate << "\n";
    std::cerr << "net_nt_mismatch_rate_per_base=" << net_nt_rate << "\n";

    // Per-codon mismatch rates: codon_mismatches / codons_compared, where a
    // codon is "compared" if at least one base in it passed its region's
    // quality filter for this read. This is a coarser view of the same
    // signal, appropriate when the DMS library mutations are codon-scale.
    double zone_codon_rate = (sum_zone_codons_compared ? (double)sum_codon_mismatches / (double)sum_zone_codons_compared : 0.0);
    double fw_codon_rate   = (sum_fw_codons_compared   ? (double)sum_fw_codon_mismatches / (double)sum_fw_codons_compared   : 0.0);
    double net_codon_rate  = zone_codon_rate - fw_codon_rate;
    if (net_codon_rate < 0) net_codon_rate = 0.0;

    std::cerr << "dms_codon_mismatch_rate_per_codon=" << zone_codon_rate << "\n";
    std::cerr << "framework_codon_mismatch_rate_per_codon=" << fw_codon_rate << "\n";
    std::cerr << "net_codon_mismatch_rate_per_codon=" << net_codon_rate << "\n";



    if (!out_counts.empty()) std::cerr << "out_counts=" << out_counts << "\n";

    // Explicit htslib cleanup: bam_destroy1 frees the record allocated by
    // bam_init1, bam_hdr_destroy frees the parsed header, sam_close flushes
    // and releases the BAM handle, fai_destroy drops the FAI index.
    bam_destroy1(b);
    bam_hdr_destroy(hdr);
    sam_close(in);
    fai_destroy(fai);
    out_fq_stream.close();

    return 0;
}
