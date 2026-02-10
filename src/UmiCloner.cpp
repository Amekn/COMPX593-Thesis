/*
Usage: UmiCloner <umi_fastq> <source_fastq> <destination_fastq>
Description:
    For fastq records in umi_fastq (with UMI in header), find matching records (by read_id) in source_fastq,
    take the header (with UMI) from umi_fastq and sequence/quality from source_fastq, and write to destination_fastq.
Author: Kang Zhou
*/

#include <iostream>
#include <fstream>
#include <string>
#include <cstdint>
#include <unordered_map>
#include <cctype>

using namespace std;

static inline char up(char c) {
    return static_cast<char>(toupper(static_cast<unsigned char>(c)));
}

// Extract read id: drop '@', take first token, then split at ':' if present
static string extract_read_id(const string& header) {
    size_t start = 0;
    size_t end = header.find_first_of(":");
    return (end == string::npos) ? string() : header.substr(start, end - start);
}

struct FastqRecord {
    string header;
    string sequence;
    string plus;
    string quality;
};

static bool read_fastq_record(istream& in, FastqRecord& r) {
    if (!getline(in, r.header)) return false;
    if (!getline(in, r.sequence)) return false;
    if (!getline(in, r.plus)) return false;
    if (!getline(in, r.quality)) return false;
    return true;
}

int main(int argc, char** argv) {
    ios::sync_with_stdio(false);
    cin.tie(nullptr);

    if (argc < 4) {
        cerr << "Usage: " << argv[0]
             << " <umi_fastq> <source_fastq> <destination_fastq>\n";
        return 1;
    }

    const string umi_fastq = argv[1];
    const string source_fastq = argv[2];
    const string destination_fastq = argv[3];

    ifstream umi_in(umi_fastq);
    if (!umi_in) {
        cerr << "ERROR: cannot open input: " << umi_fastq << "\n";
        return 2;
    }

    ifstream source_in(source_fastq);
    if (!source_in) {
        cerr << "ERROR: cannot open input: " << source_fastq << "\n";
        return 2;
    }

    ofstream out(destination_fastq);
    if (!out) {
        cerr << "ERROR: cannot open output: " << destination_fastq << "\n";
        return 2;
    }

    // read_id -> full header from umi_fastq (already contains UMI)
    unordered_map<string, string> id_to_header_map;
    id_to_header_map.reserve(1 << 24); // ~16 million entries

    int umi_total = 0, source_total = 0;
    int kept = 0, dropped = 0;

    FastqRecord record;
    record.header.reserve(128);
    record.sequence.reserve(800);
    record.plus.reserve(8);
    record.quality.reserve(800);
    string read_id;

    // Load headers from umi_fastq
    while (read_fastq_record(umi_in, record)) {
        umi_total++;

        if (record.header.empty() || record.header[0] != '@' ||
            record.plus.empty() || record.plus[0] != '+') {
            cerr << "ERROR: malformed FASTQ record in umi_fastq at record " << umi_total << "\n";
            return 3;
        }

        if (record.sequence.size() != record.quality.size()) continue;

        read_id = extract_read_id(record.header);

        if (read_id.empty()) {
            cerr << "ERROR: empty read_id in umi_fastq at record " << umi_total << "\n";
            return 3;
        }

        // Store full header
        id_to_header_map[read_id] = record.header;
    }

    // Stream through source_fastq and emit matched records
    while (read_fastq_record(source_in, record)) {
        source_total++;

        if (record.header.empty() || record.header[0] != '@' ||
            record.plus.empty() || record.plus[0] != '+') {
            cerr << "ERROR: malformed FASTQ record in source_fastq at record " << source_total << "\n";
            return 3;
        }

        if (record.sequence.size() != record.quality.size()) continue;

        auto it = id_to_header_map.find(record.header);
        if (it != id_to_header_map.end()) {
            // Use header from umi_fastq, but seq/qual from source_fastq
            out << it->second << '\n'
                << record.sequence << '\n'
                << record.plus << '\n'
                << record.quality << '\n';
            kept++;
        } else {
            dropped++;
        }
    }

    cerr << "[UmiCloner] umi_total=" << umi_total
         << " source_total=" << source_total
         << " kept=" << kept
         << " dropped=" << dropped << "\n";

    return 0;
}
