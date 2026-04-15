// VariantConcordance: population-level similarity metrics between two variant
// distributions (a "source" distribution and a "ground-truth" distribution).
//
// Library component of the thesis pipeline. Produces the static archive
// `libont_variant_concordance.a` linked by the VariantConcordance tool. Given
// two variant_key/count TSVs (as emitted by DualSiteDMSFilter --out-counts)
// and a minimum count threshold, computes: per-set haplotype counts, set
// overlap, exact overlap mass (L1 intersection of the two distributions),
// weighted Jaccard, Jensen-Shannon similarity (1 - JSD/log2), and Spearman's
// rank correlation on the top 100 most-frequent ground-truth variants.
#include "ont_tools/VariantConcordance.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <fstream>
#include <limits>
#include <numeric>
#include <stdexcept>
#include <string>
#include <unordered_set>
#include <utility>
#include <vector>

namespace ont_tools {

namespace {

// Return a copy of `text` with ASCII whitespace stripped from both ends.
// Used when parsing TSV fields where hand-edited files sometimes carry
// trailing spaces or stray CR bytes.
std::string trim_copy(const std::string& text) {
    const std::size_t start = text.find_first_not_of(" \t\r\n");
    if (start == std::string::npos) {
        return {};
    }
    const std::size_t end = text.find_last_not_of(" \t\r\n");
    return text.substr(start, end - start + 1);
}

// Split a single TSV line on literal '\t'. Preserves empty trailing fields and
// does not trim; trimming is applied per-field by the caller where needed.
std::vector<std::string> split_tsv_line(const std::string& line) {
    std::vector<std::string> fields;
    std::size_t start = 0;
    while (true) {
        const std::size_t pos = line.find('\t', start);
        if (pos == std::string::npos) {
            fields.push_back(line.substr(start));
            break;
        }
        fields.push_back(line.substr(start, pos - start));
        start = pos + 1;
    }
    return fields;
}

// Parse a variant count field as a non-negative decimal integer. Throws
// runtime_error citing `path` and `line_number` on any parse failure so the
// caller can surface a location-aware error message to the user.
std::uint64_t parse_count_or_throw(
    const std::string& text,
    const std::string& path,
    std::uint64_t line_number
) {
    const std::string trimmed = trim_copy(text);
    if (trimmed.empty()) {
        throw std::runtime_error(
            "missing count in '" + path + "' at line " + std::to_string(line_number)
        );
    }

    std::size_t consumed = 0;
    const unsigned long long value = std::stoull(trimmed, &consumed);
    if (consumed != trimmed.size()) {
        throw std::runtime_error(
            "invalid integer count '" + trimmed + "' in '" + path + "' at line "
            + std::to_string(line_number)
        );
    }
    if (value > std::numeric_limits<std::uint64_t>::max()) {
        throw std::runtime_error(
            "count out of range in '" + path + "' at line " + std::to_string(line_number)
        );
    }
    return static_cast<std::uint64_t>(value);
}

// Exact overlap mass = sum over variants of min(p_source, p_gt), i.e. the
// total probability mass both distributions agree on. This equals
// 1 - 0.5 * L1(p,q); higher is more concordant. Returns 0.0 if either input
// is empty (the distribution is undefined).
double compute_exact_overlap_mass(
    const std::unordered_map<std::string, std::uint64_t>& source_counts,
    const std::unordered_map<std::string, std::uint64_t>& ground_truth_counts
) {
    std::uint64_t source_total = 0;
    for (const auto& item : source_counts) {
        source_total += item.second;
    }

    std::uint64_t ground_truth_total = 0;
    for (const auto& item : ground_truth_counts) {
        ground_truth_total += item.second;
    }

    if (source_total == 0 || ground_truth_total == 0) {
        return 0.0;
    }

    double overlap = 0.0;
    std::unordered_set<std::string> all_sequences;
    all_sequences.reserve(source_counts.size() + ground_truth_counts.size());
    for (const auto& item : source_counts) {
        all_sequences.insert(item.first);
    }
    for (const auto& item : ground_truth_counts) {
        all_sequences.insert(item.first);
    }

    for (const auto& sequence : all_sequences) {
        const double source_prob =
            static_cast<double>(source_counts.count(sequence) ? source_counts.at(sequence) : 0ULL)
            / static_cast<double>(source_total);
        const double ground_truth_prob =
            static_cast<double>(
                ground_truth_counts.count(sequence) ? ground_truth_counts.at(sequence) : 0ULL
            ) / static_cast<double>(ground_truth_total);
        overlap += std::min(source_prob, ground_truth_prob);
    }
    return overlap;
}

// Weighted Jaccard index: sum(min(a_i, b_i)) / sum(max(a_i, b_i)) taken over
// the union of keys. Robust to the relative-abundance skew between Illumina
// and nanopore call sets, and unlike set Jaccard it penalises count
// disagreement even on shared variants. Returns 0.0 on an empty union.
double compute_weighted_jaccard(
    const std::unordered_map<std::string, std::uint64_t>& source_counts,
    const std::unordered_map<std::string, std::uint64_t>& ground_truth_counts
) {
    std::uint64_t numerator = 0;
    std::uint64_t denominator = 0;

    std::unordered_set<std::string> all_sequences;
    all_sequences.reserve(source_counts.size() + ground_truth_counts.size());
    for (const auto& item : source_counts) {
        all_sequences.insert(item.first);
    }
    for (const auto& item : ground_truth_counts) {
        all_sequences.insert(item.first);
    }

    for (const auto& sequence : all_sequences) {
        const std::uint64_t source_value =
            source_counts.count(sequence) ? source_counts.at(sequence) : 0ULL;
        const std::uint64_t ground_truth_value =
            ground_truth_counts.count(sequence) ? ground_truth_counts.at(sequence) : 0ULL;
        numerator += std::min(source_value, ground_truth_value);
        denominator += std::max(source_value, ground_truth_value);
    }

    if (denominator == 0) {
        return 0.0;
    }
    return static_cast<double>(numerator) / static_cast<double>(denominator);
}

// Jensen-Shannon similarity = 1 - JSD(p,q)/log(2), mapping a JSD in [0, log2]
// into a similarity in [0, 1] (1 == identical, 0 == disjoint). Implemented by
// summing the two KL terms against the midpoint and dividing the total by
// log(2) so the natural-log-based KL is converted to bits.
double compute_js_similarity(
    const std::unordered_map<std::string, std::uint64_t>& source_counts,
    const std::unordered_map<std::string, std::uint64_t>& ground_truth_counts
) {
    std::uint64_t source_total = 0;
    for (const auto& item : source_counts) {
        source_total += item.second;
    }

    std::uint64_t ground_truth_total = 0;
    for (const auto& item : ground_truth_counts) {
        ground_truth_total += item.second;
    }

    if (source_total == 0 || ground_truth_total == 0) {
        return 0.0;
    }

    double js_divergence = 0.0;
    std::unordered_set<std::string> all_sequences;
    all_sequences.reserve(source_counts.size() + ground_truth_counts.size());
    for (const auto& item : source_counts) {
        all_sequences.insert(item.first);
    }
    for (const auto& item : ground_truth_counts) {
        all_sequences.insert(item.first);
    }

    for (const auto& sequence : all_sequences) {
        const double source_prob =
            static_cast<double>(source_counts.count(sequence) ? source_counts.at(sequence) : 0ULL)
            / static_cast<double>(source_total);
        const double ground_truth_prob =
            static_cast<double>(
                ground_truth_counts.count(sequence) ? ground_truth_counts.at(sequence) : 0ULL
            ) / static_cast<double>(ground_truth_total);
        const double midpoint = 0.5 * (source_prob + ground_truth_prob);
        if (source_prob > 0.0) {
            js_divergence += 0.5 * source_prob * std::log(source_prob / midpoint);
        }
        if (ground_truth_prob > 0.0) {
            js_divergence += 0.5 * ground_truth_prob * std::log(ground_truth_prob / midpoint);
        }
    }

    const double normalized = js_divergence / std::log(2.0);
    return std::max(0.0, 1.0 - normalized);
}

// Compute the fractional ("average") ranks of `values`. Ties receive the
// midpoint rank of the tied block, matching the standard Spearman tie-handling
// convention. Returned vector is index-aligned with the input.
std::vector<double> average_ranks(const std::vector<double>& values) {
    std::vector<std::pair<std::size_t, double>> indexed;
    indexed.reserve(values.size());
    for (std::size_t index = 0; index < values.size(); ++index) {
        indexed.emplace_back(index, values[index]);
    }
    std::sort(indexed.begin(), indexed.end(), [](const auto& left, const auto& right) {
        return left.second < right.second;
    });

    std::vector<double> ranks(values.size(), 0.0);
    std::size_t position = 1;
    std::size_t cursor = 0;
    while (cursor < indexed.size()) {
        std::size_t end = cursor;
        const double current_value = indexed[cursor].second;
        while (end < indexed.size() && indexed[end].second == current_value) {
            ++end;
        }

        const double average_rank =
            (static_cast<double>(position)
             + static_cast<double>(position + (end - cursor) - 1))
            / 2.0;
        for (std::size_t item_index = cursor; item_index < end; ++item_index) {
            ranks[indexed[item_index].first] = average_rank;
        }

        position += end - cursor;
        cursor = end;
    }

    return ranks;
}

// Pearson correlation of two equal-length vectors. Returns 0.0 on empty/size
// mismatch. When either input has zero variance the correlation is
// conventionally undefined; we return 1.0 iff the two inputs are elementwise
// equal (within 1e-12) and 0.0 otherwise, which keeps Spearman monotonic when
// both distributions are constant.
double pearson_correlation(const std::vector<double>& left, const std::vector<double>& right) {
    if (left.size() != right.size() || left.empty()) {
        return 0.0;
    }

    const double left_mean =
        std::accumulate(left.begin(), left.end(), 0.0) / static_cast<double>(left.size());
    const double right_mean =
        std::accumulate(right.begin(), right.end(), 0.0) / static_cast<double>(right.size());

    std::vector<double> left_diff;
    left_diff.reserve(left.size());
    for (double value : left) {
        left_diff.push_back(value - left_mean);
    }

    std::vector<double> right_diff;
    right_diff.reserve(right.size());
    for (double value : right) {
        right_diff.push_back(value - right_mean);
    }

    double left_scale = 0.0;
    for (double value : left_diff) {
        left_scale += value * value;
    }
    left_scale = std::sqrt(left_scale);

    double right_scale = 0.0;
    for (double value : right_diff) {
        right_scale += value * value;
    }
    right_scale = std::sqrt(right_scale);

    if (left_scale == 0.0 || right_scale == 0.0) {
        for (std::size_t index = 0; index < left.size(); ++index) {
            if (std::fabs(left[index] - right[index]) >= 1e-12) {
                return 0.0;
            }
        }
        return 1.0;
    }

    double numerator = 0.0;
    for (std::size_t index = 0; index < left.size(); ++index) {
        numerator += left_diff[index] * right_diff[index];
    }
    return numerator / (left_scale * right_scale);
}

// Spearman rank correlation restricted to the top-100 most-abundant variants
// in the ground truth, with the corresponding source counts (0 if absent).
// Focuses the metric on the populated tail of the DMS library rather than
// letting abundant singletons dominate. Rank ties use average ranks.
double compute_top100_spearman(
    const std::unordered_map<std::string, std::uint64_t>& source_counts_all,
    const std::unordered_map<std::string, std::uint64_t>& ground_truth_high_counts
) {
    std::vector<std::pair<std::string, std::uint64_t>> ranked_ground_truth;
    ranked_ground_truth.reserve(ground_truth_high_counts.size());
    for (const auto& item : ground_truth_high_counts) {
        ranked_ground_truth.push_back(item);
    }

    // Sort by count descending; break ties by variant_key ascending so the
    // top-100 selection is deterministic across runs.
    std::sort(ranked_ground_truth.begin(), ranked_ground_truth.end(), [](const auto& left, const auto& right) {
        if (left.second != right.second) {
            return left.second > right.second;
        }
        return left.first < right.first;
    });

    if (ranked_ground_truth.empty()) {
        return 0.0;
    }

    // Single-variant GT edge case: Spearman is undefined for n=1. Treat
    // exact count equality as perfect correlation and any disagreement as
    // zero so the metric degrades gracefully rather than producing NaN.
    if (ranked_ground_truth.size() == 1U) {
        const auto& entry = ranked_ground_truth.front();
        const auto source_it = source_counts_all.find(entry.first);
        const std::uint64_t source_value = (source_it == source_counts_all.end()) ? 0ULL : source_it->second;
        return source_value == entry.second ? 1.0 : 0.0;
    }

    const std::size_t limit = std::min<std::size_t>(100, ranked_ground_truth.size());
    std::vector<double> ground_truth_values;
    std::vector<double> source_values;
    ground_truth_values.reserve(limit);
    source_values.reserve(limit);

    for (std::size_t index = 0; index < limit; ++index) {
        const auto& entry = ranked_ground_truth[index];
        ground_truth_values.push_back(static_cast<double>(entry.second));
        const auto source_it = source_counts_all.find(entry.first);
        const std::uint64_t source_value = (source_it == source_counts_all.end()) ? 0ULL : source_it->second;
        source_values.push_back(static_cast<double>(source_value));
    }

    return pearson_correlation(average_ranks(ground_truth_values), average_ranks(source_values));
}

}  // namespace

// Read a variant_key/count TSV (as emitted by DualSiteDMSFilter --out-counts)
// and aggregate duplicate keys by summing their counts. Accepts an optional
// header row (columns may be in any order provided both `variant_key` and
// `count` are present) and tolerates '#'-prefixed comment lines. Duplicate
// variant_keys in the input are additively folded into the store so that
// concatenated files or re-run outputs merge cleanly.
VariantCountStore load_variant_count_store(const std::string& path) {
    std::ifstream input(path);
    if (!input) {
        throw std::runtime_error("cannot open input file: " + path);
    }

    VariantCountStore store;
    std::string line;
    std::uint64_t line_number = 0;
    bool parsed_first_data_line = false;
    std::size_t variant_key_index = 0;
    std::size_t count_index = 1;

    while (std::getline(input, line)) {
        ++line_number;
        if (!line.empty() && line.back() == '\r') {
            line.pop_back();
        }

        const std::string trimmed = trim_copy(line);
        if (trimmed.empty() || trimmed[0] == '#') {
            continue;
        }

        const std::vector<std::string> fields = split_tsv_line(trimmed);
        if (!parsed_first_data_line) {
            parsed_first_data_line = true;
            bool header_found = false;
            for (std::size_t index = 0; index < fields.size(); ++index) {
                const std::string name = trim_copy(fields[index]);
                if (name == "variant_key") {
                    variant_key_index = index;
                    header_found = true;
                } else if (name == "count") {
                    count_index = index;
                    header_found = true;
                }
            }

            if (header_found) {
                bool has_variant_key = false;
                bool has_count = false;
                for (const auto& field : fields) {
                    const std::string name = trim_copy(field);
                    has_variant_key = has_variant_key || (name == "variant_key");
                    has_count = has_count || (name == "count");
                }
                if (!has_variant_key || !has_count) {
                    throw std::runtime_error(
                        "header in '" + path + "' must include both 'variant_key' and 'count' columns"
                    );
                }
                continue;
            }
        }

        if (fields.size() <= variant_key_index || fields.size() <= count_index) {
            throw std::runtime_error(
                "expected tab-delimited variant_key and count columns in '" + path + "' at line "
                + std::to_string(line_number)
            );
        }

        const std::string variant_key = trim_copy(fields[variant_key_index]);
        if (variant_key.empty()) {
            throw std::runtime_error(
                "empty variant_key in '" + path + "' at line " + std::to_string(line_number)
            );
        }

        const std::uint64_t count = parse_count_or_throw(fields[count_index], path, line_number);
        store.counts[variant_key] += count;
    }

    return store;
}

// Top-level entry point. Applies `haplotype_threshold` to both stores to
// filter out low-count noise (singleton errors), then computes the full
// metric battery on the filtered populations. The top-100 Spearman keeps
// operating on the unfiltered source counts so the metric is not distorted
// when a GT-abundant variant is suppressed in the source by the threshold.
ConcordanceMetrics compute_variant_concordance(
    const VariantCountStore& source_store,
    const VariantCountStore& ground_truth_store,
    const std::uint64_t haplotype_threshold
) {
    std::unordered_map<std::string, std::uint64_t> ground_truth_high_counts;
    ground_truth_high_counts.reserve(ground_truth_store.counts.size());
    for (const auto& item : ground_truth_store.counts) {
        if (item.second >= haplotype_threshold) {
            ground_truth_high_counts.insert(item);
        }
    }

    std::unordered_set<std::string> source_high_sequences;
    source_high_sequences.reserve(source_store.counts.size());
    for (const auto& item : source_store.counts) {
        if (item.second >= haplotype_threshold) {
            source_high_sequences.insert(item.first);
        }
    }

    std::unordered_set<std::string> union_sequences = source_high_sequences;
    union_sequences.reserve(source_high_sequences.size() + ground_truth_high_counts.size());
    for (const auto& item : ground_truth_high_counts) {
        union_sequences.insert(item.first);
    }

    std::size_t intersection_haplotype_count = 0;
    for (const auto& item : ground_truth_high_counts) {
        if (source_high_sequences.find(item.first) != source_high_sequences.end()) {
            ++intersection_haplotype_count;
        }
    }

    // Build dense count vectors over the union of threshold-passing sequences.
    // For each sequence, source-only or gt-only keys receive a zero count on
    // the opposite side, which is exactly the semantics the mass/JSD/weighted
    // Jaccard math requires.
    std::unordered_map<std::string, std::uint64_t> filtered_source_counts;
    filtered_source_counts.reserve(union_sequences.size());
    std::unordered_map<std::string, std::uint64_t> filtered_ground_truth_counts;
    filtered_ground_truth_counts.reserve(union_sequences.size());

    for (const auto& sequence : union_sequences) {
        const auto source_it = source_store.counts.find(sequence);
        filtered_source_counts[sequence] =
            (source_it == source_store.counts.end()) ? 0ULL : source_it->second;

        const auto ground_truth_it = ground_truth_store.counts.find(sequence);
        filtered_ground_truth_counts[sequence] =
            (ground_truth_it == ground_truth_store.counts.end()) ? 0ULL : ground_truth_it->second;
    }

    ConcordanceMetrics metrics;
    metrics.haplotype_threshold = haplotype_threshold;
    metrics.source_haplotypes = source_high_sequences.size();
    metrics.groundtruth_haplotypes = ground_truth_high_counts.size();
    metrics.intersection_haplotypes = intersection_haplotype_count;
    metrics.union_haplotypes = union_sequences.size();
    metrics.exact_overlap_mass =
        compute_exact_overlap_mass(filtered_source_counts, filtered_ground_truth_counts);
    metrics.weighted_jaccard =
        compute_weighted_jaccard(filtered_source_counts, filtered_ground_truth_counts);
    metrics.jensen_shannon_similarity =
        compute_js_similarity(filtered_source_counts, filtered_ground_truth_counts);
    metrics.top100_spearman =
        compute_top100_spearman(source_store.counts, ground_truth_high_counts);
    return metrics;
}

}  // namespace ont_tools
