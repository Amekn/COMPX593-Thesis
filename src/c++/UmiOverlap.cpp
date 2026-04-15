// UmiOverlap.cpp
// Produces executable: UmiOverlap
// Usage: UmiOverlap <fastq1> <fastq2>
//
// Pipeline role: compute unique-UMI-set overlap between two FASTQ files
// annotated with ':UMI_<key>' tags. Used to gauge how much the nanopore
// and MGI-derived UMI populations share after dedup, and to sanity-check
// that UmiFilter output has the expected cardinality.
#include <cstdint>
#include <fstream>
#include <iostream>
#include <string>
#include <unordered_set>

#include "ont_tools/Fastq.hpp"
#include "ont_tools/Umi.hpp"

namespace {

// Print CLI contract to stderr when argv count does not match.
void PrintUsage(const char* program_name) {
    std::cerr
        << "Usage: " << program_name << " <fastq1> <fastq2>\n\n"
        << "Report read counts and the overlap between the UMI key sets of two FASTQ files.\n";
}

// Per-file summary: total records plus the set of distinct UMI keys seen.
struct FileSummary {
    std::uint64_t read_count = 0;
    std::unordered_set<std::string> umi_keys;
};

// Stream a FASTQ and collect distinct UMI keys plus the record count.
// Throws std::runtime_error if the path cannot be opened.
FileSummary LoadUmiKeys(const std::string& fastq_path) {
    FileSummary summary;
    // ~4M buckets: matches the UmiFilter sizing so a deduplicated FASTQ
    // can be consumed without rehashing on the expected-case library size.
    summary.umi_keys.reserve(1U << 22);

    std::ifstream input_stream(fastq_path);
    if (!input_stream.is_open()) {
        throw std::runtime_error("Failed to open FASTQ: " + fastq_path);
    }

    ont::fastq::Record record;
    while (ont::fastq::ReadRecord(input_stream, record, fastq_path)) {
        ++summary.read_count;
        summary.umi_keys.insert(ont::umi::ParseSingleKey(record.header));
    }

    return summary;
}

// Count how many keys are present in both sets. Iterates the smaller set and
// probes the larger, so cost is O(min(|a|,|b|)) average-case.
std::uint64_t CountSetOverlap(
    const std::unordered_set<std::string>& left_keys,
    const std::unordered_set<std::string>& right_keys) {
    const auto* smaller_set = &left_keys;
    const auto* larger_set = &right_keys;
    if (left_keys.size() > right_keys.size()) {
        smaller_set = &right_keys;
        larger_set = &left_keys;
    }

    std::uint64_t overlap_count = 0;
    for (const std::string& key : *smaller_set) {
        if (larger_set->find(key) != larger_set->end()) {
            ++overlap_count;
        }
    }
    return overlap_count;
}

}  // namespace

int main(int argc, char** argv) {
    std::ios::sync_with_stdio(false);
    std::cin.tie(nullptr);

    if (argc != 3) {
        PrintUsage(argv[0]);
        return 1;
    }

    try {
        const FileSummary first_summary = LoadUmiKeys(argv[1]);
        const FileSummary second_summary = LoadUmiKeys(argv[2]);
        const std::uint64_t overlapping_unique_umis =
            CountSetOverlap(first_summary.umi_keys, second_summary.umi_keys);

        std::cout << "fastq1_reads=" << first_summary.read_count
                  << " fastq2_reads=" << second_summary.read_count
                  << " fastq1_unique_umis=" << first_summary.umi_keys.size()
                  << " fastq2_unique_umis=" << second_summary.umi_keys.size()
                  << " overlapping_unique_umis=" << overlapping_unique_umis
                  << "\n";
    } catch (const std::exception& exception) {
        std::cerr << "Error: " << exception.what() << "\n";
        return 2;
    }

    return 0;
}
