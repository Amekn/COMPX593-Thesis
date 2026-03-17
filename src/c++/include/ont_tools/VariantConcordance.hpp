#pragma once

#include <cstddef>
#include <cstdint>
#include <string>
#include <unordered_map>

namespace ont_tools {

struct VariantCountStore {
    std::unordered_map<std::string, std::uint64_t> counts;
};

struct ConcordanceMetrics {
    std::uint64_t haplotype_threshold = 0;
    std::size_t source_haplotypes = 0;
    std::size_t groundtruth_haplotypes = 0;
    std::size_t intersection_haplotypes = 0;
    std::size_t union_haplotypes = 0;
    double exact_overlap_mass = 0.0;
    double weighted_jaccard = 0.0;
    double jensen_shannon_similarity = 0.0;
    double top100_spearman = 0.0;
};

VariantCountStore load_variant_count_store(const std::string& path);

ConcordanceMetrics compute_variant_concordance(
    const VariantCountStore& source_store,
    const VariantCountStore& ground_truth_store,
    std::uint64_t haplotype_threshold
);

}  // namespace ont_tools
