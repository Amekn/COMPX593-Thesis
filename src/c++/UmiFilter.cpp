#include <cstdint>
#include <iostream>
#include <string>
#include <unordered_set>

#include "ont_tools/Fastq.hpp"
#include "ont_tools/Umi.hpp"

namespace {

void PrintUsage(const char* program_name) {
    std::cerr
        << "Usage: " << program_name << " <input.fastq> <output.fastq>\n\n"
        << "Keep only the first read observed for each UMI key in the FASTQ header.\n";
}

}  // namespace

int main(int argc, char** argv) {
    std::ios::sync_with_stdio(false);
    std::cin.tie(nullptr);

    if (argc != 3) {
        PrintUsage(argv[0]);
        return 1;
    }

    const std::string input_fastq = argv[1];
    const std::string output_fastq = argv[2];

    std::unordered_set<std::string> seen_umi_keys;
    seen_umi_keys.reserve(1U << 22);

    std::uint64_t total_reads = 0;
    std::uint64_t kept_reads = 0;
    std::uint64_t duplicate_umi_reads = 0;
    std::uint64_t missing_umi_reads = 0;

    try {
        ont::fastq::Reader reader(input_fastq);
        ont::fastq::Writer writer(output_fastq);

        ont::fastq::Record record;
        while (reader.Next(record)) {
            ++total_reads;

            std::string umi_key;
            try {
                umi_key = ont::umi::ParseSingleKey(record.header);
            } catch (const std::exception&) {
                ++missing_umi_reads;
                continue;
            }

            const auto [_, inserted] = seen_umi_keys.insert(std::move(umi_key));
            if (inserted) {
                writer.Write(record);
                ++kept_reads;
            } else {
                ++duplicate_umi_reads;
            }
        }
    } catch (const std::exception& exception) {
        std::cerr << "Error: " << exception.what() << "\n";
        return 2;
    }

    std::cerr << "[UmiFilter] total=" << total_reads
              << " kept=" << kept_reads
              << " duplicate_umi_reads=" << duplicate_umi_reads
              << " missing_umi_reads=" << missing_umi_reads << "\n";
    return 0;
}
