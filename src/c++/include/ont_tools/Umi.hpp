// Umi.hpp
// Helpers for parsing UMI tags out of FASTQ headers produced by the nested
// PCR protocol. The pipeline appends ":UMI_<fwd>_<rev>" to the read name so
// nanopore and MGI reads can be cross-referenced. Routines here extract the
// read name, the whole tag value, or the split dual-UMI halves; all parsing
// is header-only so each tool compiles it directly without a library.
#pragma once

#include <stdexcept>
#include <string>
#include <string_view>

#include "ont_tools/Sequence.hpp"

namespace ont::umi {

// Split representation of a dual UMI parsed from a header tag.
struct DualKey {
    std::string first;
    std::string second;
};

// Return the first whitespace-delimited token of a FASTQ header, with the
// leading '@' stripped if present. Accepts either a raw or '@'-prefixed view.
inline std::string_view HeaderToken(std::string_view header) {
    if (!header.empty() && header.front() == '@') {
        header.remove_prefix(1);
    }

    const std::size_t token_end = header.find_first_of(" \t\r\n");
    return header.substr(0, token_end);
}

// Return the read name with any ":UMI_..." tag stripped. Used to match the
// same read across UMI-annotated and vanilla FASTQ streams (e.g. in
// UmiCloner) where the tag is the only difference between headers.
inline std::string ExtractReadName(std::string_view header) {
    std::string_view token = HeaderToken(header);
    const std::size_t umi_tag_position = token.find(":UMI_");
    if (umi_tag_position != std::string_view::npos) {
        token = token.substr(0, umi_tag_position);
    }
    return std::string(token);
}

// Return the upper-cased ":UMI_" value (everything between the tag and the
// next whitespace). Throws std::runtime_error if the tag is absent or empty.
inline std::string ExtractUmiTagValue(std::string_view header) {
    constexpr std::string_view kUmiTag = ":UMI_";
    const std::size_t tag_position = header.find(kUmiTag);
    if (tag_position == std::string_view::npos) {
        throw std::runtime_error("Header is missing the ':UMI_' tag.");
    }

    std::size_t value_begin = tag_position + kUmiTag.size();
    std::size_t value_end = header.find_first_of(" \t\r\n", value_begin);
    if (value_end == std::string_view::npos) {
        value_end = header.size();
    }
    if (value_end <= value_begin) {
        throw std::runtime_error("Header contains ':UMI_' without a value.");
    }

    return ont::seq::ToUpperCopy(std::string(header.substr(value_begin, value_end - value_begin)));
}

// Return the whole UMI tag value as a single string; used by tools that
// treat a dual UMI as an atomic key (e.g. UmiFilter single-pass dedup).
inline std::string ParseSingleKey(std::string_view header) {
    return ExtractUmiTagValue(header);
}

// Split the UMI tag value at the first underscore into its forward and
// reverse halves. Throws std::runtime_error if the separator is missing or
// either side is empty.
inline DualKey ParseDualKey(std::string_view header) {
    const std::string value = ExtractUmiTagValue(header);
    const std::size_t separator = value.find('_');
    if (separator == std::string::npos) {
        throw std::runtime_error("Dual UMI header is missing the underscore separator.");
    }

    DualKey key{value.substr(0, separator), value.substr(separator + 1)};
    if (key.first.empty() || key.second.empty()) {
        throw std::runtime_error("Dual UMI header contains an empty UMI component.");
    }
    return key;
}

// Parse and re-join the dual UMI with an underscore, yielding a canonical
// "fwd_rev" hash key independent of any surrounding header whitespace.
inline std::string ParseDualKeyString(std::string_view header) {
    const DualKey key = ParseDualKey(header);
    return key.first + "_" + key.second;
}

}  // namespace ont::umi
