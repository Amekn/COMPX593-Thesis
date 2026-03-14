#pragma once

#include <stdexcept>
#include <string>
#include <string_view>

#include "ont_tools/Sequence.hpp"

namespace ont::umi {

struct DualKey {
    std::string first;
    std::string second;
};

inline std::string_view HeaderToken(std::string_view header) {
    if (!header.empty() && header.front() == '@') {
        header.remove_prefix(1);
    }

    const std::size_t token_end = header.find_first_of(" \t\r\n");
    return header.substr(0, token_end);
}

inline std::string ExtractReadName(std::string_view header) {
    std::string_view token = HeaderToken(header);
    const std::size_t umi_tag_position = token.find(":UMI_");
    if (umi_tag_position != std::string_view::npos) {
        token = token.substr(0, umi_tag_position);
    }
    return std::string(token);
}

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

inline std::string ParseSingleKey(std::string_view header) {
    return ExtractUmiTagValue(header);
}

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

inline std::string ParseDualKeyString(std::string_view header) {
    const DualKey key = ParseDualKey(header);
    return key.first + "_" + key.second;
}

}  // namespace ont::umi
