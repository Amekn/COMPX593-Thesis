#include <cstdint>
#include <iostream>
#include <string>
#include <unordered_map>

#include "ont_tools/Fastq.hpp"
#include "ont_tools/Umi.hpp"

namespace {

void PrintUsage(const char* program_name) {
    std::cerr
        << "Usage: " << program_name << " <umi_fastq> <source_fastq> <destination_fastq>\n\n"
        << "Copy UMI-enriched headers from the first FASTQ onto matching reads from the second FASTQ.\n"
        << "Read matching is performed on the FASTQ read name, independent of any appended UMI tag.\n";
}

}  // namespace

int main(int argc, char** argv) {
    std::ios::sync_with_stdio(false);
    std::cin.tie(nullptr);

    if (argc != 4) {
        PrintUsage(argv[0]);
        return 1;
    }

    const std::string umi_fastq = argv[1];
    const std::string source_fastq = argv[2];
    const std::string destination_fastq = argv[3];

    std::unordered_map<std::string, std::string> umi_header_by_read_name;
    umi_header_by_read_name.reserve(1U << 24);

    std::uint64_t umi_records = 0;
    std::uint64_t duplicate_umi_headers = 0;
    std::uint64_t source_records = 0;
    std::uint64_t cloned_records = 0;
    std::uint64_t unmatched_records = 0;

    try {
        ont::fastq::Reader umi_reader(umi_fastq);
        ont::fastq::Record record;
        while (umi_reader.Next(record)) {
            ++umi_records;
            const std::string read_name = ont::umi::ExtractReadName(record.header);
            const auto [_, inserted] = umi_header_by_read_name.emplace(read_name, record.header);
            if (!inserted) {
                ++duplicate_umi_headers;
            }
        }

        ont::fastq::Reader source_reader(source_fastq);
        ont::fastq::Writer destination_writer(destination_fastq);

        ont::fastq::Record source_record;
        while (source_reader.Next(source_record)) {
            ++source_records;
            const std::string read_name = ont::umi::ExtractReadName(source_record.header);
            const auto match = umi_header_by_read_name.find(read_name);
            if (match == umi_header_by_read_name.end()) {
                ++unmatched_records;
                continue;
            }

            ont::fastq::Record output_record = source_record;
            output_record.header = match->second;
            destination_writer.Write(output_record);
            ++cloned_records;
        }
    } catch (const std::exception& exception) {
        std::cerr << "Error: " << exception.what() << "\n";
        return 2;
    }

    std::cerr << "[UmiCloner] umi_records=" << umi_records
              << " duplicate_umi_headers=" << duplicate_umi_headers
              << " source_records=" << source_records
              << " cloned_records=" << cloned_records
              << " unmatched_records=" << unmatched_records << "\n";
    return 0;
}
