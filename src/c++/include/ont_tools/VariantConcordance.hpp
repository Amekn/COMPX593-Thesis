// VariantConcordance.hpp
// Declares the variant-count store and concordance-metrics data shared by
// the nanopore vs short-read haplotype comparison. Exposes a loader for
// variant_key/count TSV files and a single entry point that computes the
// thresholded concordance summary consumed by VariantConcordanceCli.
#pragma once

#include <cstddef>
#include <cstdint>
#include <string>
#include <unordered_map>

namespace ont_tools {

// Flat variant_key -> observation count mapping, one store per input TSV.
struct VariantCountStore {
    std::unordered_map<std::string, std::uint64_t> counts;
};

// Aggregate concordance figures reported for a single threshold pass; values
// are computed jointly so downstream scripts can emit a single-row TSV.
struct ConcordanceMetrics {
    std::uint64_t haplotype_threshold = 0;      // Minimum count for a haplotype to be considered "observed".
    std::size_t source_haplotypes = 0;          // Nanopore haplotypes with count >= threshold.
    std::size_t groundtruth_haplotypes = 0;     // Ground-truth haplotypes with count >= threshold.
    std::size_t intersection_haplotypes = 0;    // Haplotypes above threshold in both stores.
    std::size_t union_haplotypes = 0;           // Haplotypes above threshold in either store.
    double exact_overlap_mass = 0.0;            // Fraction of source mass whose keys also appear in ground truth.
    double weighted_jaccard = 0.0;              // Count-weighted Jaccard similarity across the union.
    double jensen_shannon_similarity = 0.0;     // 1 - sqrt(JSD) on normalised count distributions.
    double top100_spearman = 0.0;               // Spearman rank correlation on the ground-truth top-100 haplotypes.
};

// Load a TSV with "variant_key<TAB>count" rows into a VariantCountStore.
// Throws std::runtime_error on unreadable files or malformed rows.
VariantCountStore load_variant_count_store(const std::string& path);

// Compute all ConcordanceMetrics fields in one pass over the two stores,
// applying haplotype_threshold to both. The metric definitions match
// experimental/ont/scripts/model_performance_report.py.
ConcordanceMetrics compute_variant_concordance(
    const VariantCountStore& source_store,
    const VariantCountStore& ground_truth_store,
    std::uint64_t haplotype_threshold
);

}  // namespace ont_tools
