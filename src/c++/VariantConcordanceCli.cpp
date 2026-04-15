// VariantConcordanceCli.cpp
// Produces executable: VariantConcordance
// Usage: VariantConcordance <source_variant.tsv> <ground_truth_variant.tsv> <haplotype_threshold>
//
// Pipeline role: top-level driver that computes exact-haplotype concordance
// between a nanopore variant_key/count TSV and a ground-truth TSV (typically
// MGI short-read). Loads both stores via ont_tools and writes a one-row
// TSV with the metrics defined in ConcordanceMetrics. Metric behaviour is
// aligned with experimental/ont/scripts/model_performance_report.py.
#include "ont_tools/VariantConcordance.hpp"

#include <cstdint>
#include <exception>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <string>

using namespace std;

// Print CLI contract and metric documentation to stderr.
static void usage(const char* prog) {
    cerr << "Usage: " << prog << " <source_variant.tsv> <ground_truth_variant.tsv> <haplotype_threshold>\n\n"
         <<
R"(Description:
  Compute exact-haplotype concordance metrics from two variant_key/count TSV stores.

Input format:
  Expected TSV columns: variant_key<TAB>count

Output columns:
  haplotype_threshold
  source_haplotypes
    Source haplotypes with count >= haplotype_threshold.
  groundtruth_haplotypes
    Ground-truth haplotypes with count >= haplotype_threshold.
  intersection_haplotypes
    Shared haplotypes above threshold in both source and ground truth.
  union_haplotypes
    Union of source-above-threshold haplotypes and ground-truth-above-threshold haplotypes.
  exact_overlap_mass
  weighted_jaccard
  jensen_shannon_similarity
  top100_spearman

Metric behavior follows experimental/ont/scripts/model_performance_report.py
for count-based concordance, with top100_spearman ranked from the ground-truth
top 100 haplotypes above threshold.
)";
}

// Parse the CLI threshold argument as an unsigned 64-bit integer. Rejects
// strings with trailing garbage so e.g. "5x" does not silently pass as 5.
static std::uint64_t parse_threshold(const string& text) {
    size_t consumed = 0;
    const unsigned long long value = stoull(text, &consumed);
    if (consumed != text.size()) {
        throw runtime_error("invalid haplotype_threshold: " + text);
    }
    return static_cast<std::uint64_t>(value);
}

int main(int argc, char** argv) {
    ios::sync_with_stdio(false);
    cin.tie(nullptr);

    try {
        if (argc == 2) {
            const string arg = argv[1];
            if (arg == "--help" || arg == "-h") {
                usage(argv[0]);
                return 0;
            }
        }

        if (argc != 4) {
            usage(argv[0]);
            return 1;
        }

        const string source_path = argv[1];
        const string ground_truth_path = argv[2];
        const std::uint64_t haplotype_threshold = parse_threshold(argv[3]);

        // Load both stores into memory before metric computation: the
        // metrics require global lookups on both sides.
        const auto source_store = ont_tools::load_variant_count_store(source_path);
        const auto ground_truth_store = ont_tools::load_variant_count_store(ground_truth_path);
        const auto metrics = ont_tools::compute_variant_concordance(
            source_store,
            ground_truth_store,
            haplotype_threshold
        );

        // Emit a two-row TSV: a header followed by the single metric row.
        // setprecision(10) is used for fractional fields to keep downstream
        // R/Python parsers from silently rounding.
        cout << "haplotype_threshold\tsource_haplotypes\tgroundtruth_haplotypes\t"
             << "intersection_haplotypes\tunion_haplotypes\texact_overlap_mass\t"
             << "weighted_jaccard\tjensen_shannon_similarity\ttop100_spearman\n";

        cout << fixed << setprecision(10)
             << metrics.haplotype_threshold << '\t'
             << metrics.source_haplotypes << '\t'
             << metrics.groundtruth_haplotypes << '\t'
             << metrics.intersection_haplotypes << '\t'
             << metrics.union_haplotypes << '\t'
             << metrics.exact_overlap_mass << '\t'
             << metrics.weighted_jaccard << '\t'
             << metrics.jensen_shannon_similarity << '\t'
             << metrics.top100_spearman << '\n';

        return 0;
    } catch (const exception& error) {
        cerr << "Error: " << error.what() << '\n';
        return 2;
    }
}
