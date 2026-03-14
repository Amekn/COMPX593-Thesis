#include <cstdint>
#include <fstream>
#include <iostream>
#include <string>
#include <unordered_set>

#include "ont_tools/Fastq.hpp"
#include "ont_tools/Umi.hpp"

namespace {

void PrintUsage(const char* program_name) {
    std::cerr
        << "Usage: " << program_name << " <fastq1> <fastq2>\n\n"
        << "Report read counts and the overlap between the UMI key sets of two FASTQ files.\n";
}

struct FileSummary {
    std::uint64_t read_count = 0;
    std::unordered_set<std::string> umi_keys;
};

FileSummary LoadUmiKeys(const std::string& fastq_path) {
    FileSummary summary;
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
