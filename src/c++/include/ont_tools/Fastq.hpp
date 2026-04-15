// Fastq.hpp
// Minimal header-only FASTQ record type plus plain-stream Reader/Writer
// helpers. The thesis pipeline is happy with uncompressed text FASTQ
// (compression is handled externally by pigz/zcat), so the implementation
// avoids any library dependency and delegates gzip/bgzf concerns to the
// shell glue.
#pragma once

#include <fstream>
#include <istream>
#include <ostream>
#include <stdexcept>
#include <string>

namespace ont::fastq {

// A single FASTQ record held as four line-strings; plus_line is kept verbatim
// so round-trips preserve any tag appended after the '+'.
struct Record {
    std::string header;
    std::string sequence;
    std::string plus_line;
    std::string quality;
};

// Strip a trailing '\r' so Windows-line-ending FASTQs parse cleanly.
inline void TrimTrailingCarriageReturn(std::string* text) {
    if (!text->empty() && text->back() == '\r') {
        text->pop_back();
    }
}

// Read the next four-line FASTQ record from input into record. Returns false
// at EOF on the first line. Throws std::runtime_error on any structural
// failure (missing line, bad prefix, sequence/quality length mismatch);
// source_name is embedded in the exception for diagnostic context.
inline bool ReadRecord(std::istream& input, Record& record, const std::string& source_name = "FASTQ") {
    if (!std::getline(input, record.header)) {
        return false;
    }
    if (!std::getline(input, record.sequence)) {
        throw std::runtime_error("Malformed " + source_name + ": missing sequence line.");
    }
    if (!std::getline(input, record.plus_line)) {
        throw std::runtime_error("Malformed " + source_name + ": missing '+' line.");
    }
    if (!std::getline(input, record.quality)) {
        throw std::runtime_error("Malformed " + source_name + ": missing quality line.");
    }

    TrimTrailingCarriageReturn(&record.header);
    TrimTrailingCarriageReturn(&record.sequence);
    TrimTrailingCarriageReturn(&record.plus_line);
    TrimTrailingCarriageReturn(&record.quality);

    if (record.header.empty() || record.header.front() != '@') {
        throw std::runtime_error("Malformed " + source_name + ": header does not start with '@'.");
    }
    if (record.plus_line.empty() || record.plus_line.front() != '+') {
        throw std::runtime_error("Malformed " + source_name + ": third line does not start with '+'.");
    }
    if (record.sequence.size() != record.quality.size()) {
        throw std::runtime_error(
            "Malformed " + source_name + ": sequence and quality lengths differ (" +
            std::to_string(record.sequence.size()) + " vs " + std::to_string(record.quality.size()) + ").");
    }

    return true;
}

// Owning FASTQ input handle. Constructor opens the file (throws on failure)
// and retains the path so ReadRecord can surface it in error messages.
class Reader {
public:
    explicit Reader(const std::string& path) : path_(path), stream_(path) {
        if (!stream_.is_open()) {
            throw std::runtime_error("Failed to open input FASTQ: " + path);
        }
    }

    bool Next(Record& record) {
        return ReadRecord(stream_, record, path_);
    }

private:
    std::string path_;
    std::ifstream stream_;
};

// Owning FASTQ output handle. Writes records in canonical four-line form
// with a trailing '\n' on each line; no buffering tweaks beyond std::ofstream.
class Writer {
public:
    explicit Writer(const std::string& path) : stream_(path) {
        if (!stream_.is_open()) {
            throw std::runtime_error("Failed to open output FASTQ: " + path);
        }
    }

    void Write(const Record& record) {
        stream_ << record.header << '\n'
                << record.sequence << '\n'
                << record.plus_line << '\n'
                << record.quality << '\n';
    }

private:
    std::ofstream stream_;
};

}  // namespace ont::fastq
