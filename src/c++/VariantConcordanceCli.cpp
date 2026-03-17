#include "ont_tools/VariantConcordance.hpp"

#include <cstdint>
#include <exception>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <string>

using namespace std;

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

        const auto source_store = ont_tools::load_variant_count_store(source_path);
        const auto ground_truth_store = ont_tools::load_variant_count_store(ground_truth_path);
        const auto metrics = ont_tools::compute_variant_concordance(
            source_store,
            ground_truth_store,
            haplotype_threshold
        );

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
