#include <cstdint>
#include <iostream>
#include <string>
#include <unordered_set>

#include "ont_tools/Fastq.hpp"
#include "ont_tools/Umi.hpp"

namespace {

void PrintUsage(const char* program_name) {
    std::cerr
        << "Usage: " << program_name << " <fastq1> <fastq2>\n\n"
        << "Summarise read-level and unique-key overlap between two single-UMI FASTQ files.\n";
}

struct OverlapSummary {
    std::uint64_t read_count = 0;
    std::unordered_set<std::string> umi_keys;
};

OverlapSummary LoadSingleUmiFile(const std::string& fastq_path) {
    OverlapSummary summary;
    summary.umi_keys.reserve(3'000'000);
    summary.umi_keys.max_load_factor(0.70F);

    ont::fastq::Reader reader(fastq_path);
    ont::fastq::Record record;
    while (reader.Next(record)) {
        ++summary.read_count;
        summary.umi_keys.insert(ont::umi::ParseSingleKey(record.header));
    }

    return summary;
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
        const OverlapSummary first_summary = LoadSingleUmiFile(argv[1]);
        const OverlapSummary second_summary = LoadSingleUmiFile(argv[2]);

        std::unordered_set<std::string> seen_in_second_file;
        seen_in_second_file.reserve(3'000'000);
        seen_in_second_file.max_load_factor(0.70F);

        ont::fastq::Reader second_reader(argv[2]);
        ont::fastq::Record record;

        std::uint64_t overlapping_reads_in_second_file = 0;
        std::uint64_t second_only_unique_umis = 0;
        std::uint64_t overlapping_unique_umis = 0;

        while (second_reader.Next(record)) {
            const std::string umi_key = ont::umi::ParseSingleKey(record.header);
            const bool present_in_first_file =
                first_summary.umi_keys.find(umi_key) != first_summary.umi_keys.end();
            if (present_in_first_file) {
                ++overlapping_reads_in_second_file;
            }

            if (seen_in_second_file.insert(umi_key).second) {
                if (present_in_first_file) {
                    ++overlapping_unique_umis;
                } else {
                    ++second_only_unique_umis;
                }
            }
        }

        const std::uint64_t first_only_unique_umis =
            first_summary.umi_keys.size() >= overlapping_unique_umis
                ? first_summary.umi_keys.size() - overlapping_unique_umis
                : 0;

        std::cout << "FASTQ1_Reads=" << first_summary.read_count
                  << " FASTQ2_Reads=" << second_summary.read_count
                  << " FASTQ1_UniqueUMIs=" << first_summary.umi_keys.size()
                  << " FASTQ1_UniqueUMIsOnly=" << first_only_unique_umis
                  << " FASTQ2_UniqueUMIsOnly=" << second_only_unique_umis
                  << " Overlapping_UniqueUMIs=" << overlapping_unique_umis
                  << " FASTQ2_OverlapReads=" << overlapping_reads_in_second_file
                  << "\n";
    } catch (const std::exception& exception) {
        std::cerr << "Error: " << exception.what() << "\n";
        return 2;
    }

    return 0;
}
