#pragma once

#include <fstream>
#include <istream>
#include <ostream>
#include <stdexcept>
#include <string>

namespace ont::fastq {

struct Record {
    std::string header;
    std::string sequence;
    std::string plus_line;
    std::string quality;
};

inline void TrimTrailingCarriageReturn(std::string* text) {
    if (!text->empty() && text->back() == '\r') {
        text->pop_back();
    }
}

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
