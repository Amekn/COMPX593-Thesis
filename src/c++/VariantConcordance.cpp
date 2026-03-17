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

std::string trim_copy(const std::string& text) {
    const std::size_t start = text.find_first_not_of(" \t\r\n");
    if (start == std::string::npos) {
        return {};
    }
    const std::size_t end = text.find_last_not_of(" \t\r\n");
    return text.substr(start, end - start + 1);
}

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

double compute_top100_spearman(
    const std::unordered_map<std::string, std::uint64_t>& source_counts_all,
    const std::unordered_map<std::string, std::uint64_t>& ground_truth_high_counts
) {
    std::vector<std::pair<std::string, std::uint64_t>> ranked_ground_truth;
    ranked_ground_truth.reserve(ground_truth_high_counts.size());
    for (const auto& item : ground_truth_high_counts) {
        ranked_ground_truth.push_back(item);
    }

    std::sort(ranked_ground_truth.begin(), ranked_ground_truth.end(), [](const auto& left, const auto& right) {
        if (left.second != right.second) {
            return left.second > right.second;
        }
        return left.first < right.first;
    });

    if (ranked_ground_truth.empty()) {
        return 0.0;
    }

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
