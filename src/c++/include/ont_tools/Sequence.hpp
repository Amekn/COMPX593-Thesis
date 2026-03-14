#pragma once

#include <algorithm>
#include <array>
#include <cctype>
#include <cstdint>
#include <string>
#include <string_view>
#include <vector>

namespace ont::seq {

inline char ToUpper(char value) {
    return static_cast<char>(std::toupper(static_cast<unsigned char>(value)));
}

inline bool IsCanonicalBase(char value) {
    const char base = ToUpper(value);
    return base == 'A' || base == 'C' || base == 'G' || base == 'T';
}

namespace detail {

inline const std::array<std::uint8_t, 256>& BaseMaskTable() {
    static const std::array<std::uint8_t, 256> table = [] {
        std::array<std::uint8_t, 256> masks{};

        auto set_mask = [&masks](char code, const char* allowed_bases) {
            std::uint8_t mask = 0;
            for (const char* base = allowed_bases; *base != '\0'; ++base) {
                switch (*base) {
                    case 'A':
                        mask |= 1;
                        break;
                    case 'C':
                        mask |= 2;
                        break;
                    case 'G':
                        mask |= 4;
                        break;
                    case 'T':
                        mask |= 8;
                        break;
                    default:
                        break;
                }
            }
            masks[static_cast<unsigned char>(code)] = mask;
        };

        set_mask('A', "A");
        set_mask('C', "C");
        set_mask('G', "G");
        set_mask('T', "T");
        set_mask('R', "AG");
        set_mask('Y', "CT");
        set_mask('S', "GC");
        set_mask('W', "AT");
        set_mask('K', "GT");
        set_mask('M', "AC");
        set_mask('B', "CGT");
        set_mask('D', "AGT");
        set_mask('H', "ACT");
        set_mask('V', "ACG");
        set_mask('N', "ACGT");

        for (int character_code = 0; character_code < 256; ++character_code) {
            const char uppercase_code = ToUpper(static_cast<char>(character_code));
            masks[static_cast<unsigned char>(character_code)] =
                masks[static_cast<unsigned char>(uppercase_code)];
        }

        return masks;
    }();

    return table;
}

inline const std::array<char, 256>& ComplementTable() {
    static const std::array<char, 256> table = [] {
        std::array<char, 256> complements{};
        complements.fill('N');

        auto set_complement = [&complements](char input, char output) {
            complements[static_cast<unsigned char>(input)] = output;
        };

        set_complement('A', 'T');
        set_complement('C', 'G');
        set_complement('G', 'C');
        set_complement('T', 'A');
        set_complement('R', 'Y');
        set_complement('Y', 'R');
        set_complement('S', 'S');
        set_complement('W', 'W');
        set_complement('K', 'M');
        set_complement('M', 'K');
        set_complement('B', 'V');
        set_complement('D', 'H');
        set_complement('H', 'D');
        set_complement('V', 'B');
        set_complement('N', 'N');

        for (int character_code = 0; character_code < 256; ++character_code) {
            const char uppercase_code = ToUpper(static_cast<char>(character_code));
            complements[static_cast<unsigned char>(character_code)] =
                complements[static_cast<unsigned char>(uppercase_code)];
        }

        return complements;
    }();

    return table;
}

}  // namespace detail

inline std::uint8_t BaseMask(char value) {
    return detail::BaseMaskTable()[static_cast<unsigned char>(value)];
}

inline char Complement(char value) {
    return detail::ComplementTable()[static_cast<unsigned char>(value)];
}

inline std::string ToUpperCopy(std::string value) {
    std::transform(value.begin(), value.end(), value.begin(), ToUpper);
    return value;
}

inline void UppercaseInPlace(std::string* value) {
    std::transform(value->begin(), value->end(), value->begin(), ToUpper);
}

inline std::string ReverseCopy(std::string_view value) {
    return std::string(value.rbegin(), value.rend());
}

inline std::string ReverseQuality(std::string_view quality) {
    return ReverseCopy(quality);
}

inline std::string ReverseComplement(std::string_view sequence) {
    std::string reverse_complement;
    reverse_complement.reserve(sequence.size());
    for (auto iterator = sequence.rbegin(); iterator != sequence.rend(); ++iterator) {
        reverse_complement.push_back(Complement(*iterator));
    }
    return reverse_complement;
}

inline std::vector<int> AmbiguousPositions(std::string_view primer) {
    std::vector<int> positions;
    positions.reserve(primer.size());
    for (std::size_t index = 0; index < primer.size(); ++index) {
        if (!IsCanonicalBase(primer[index])) {
            positions.push_back(static_cast<int>(index));
        }
    }
    return positions;
}

inline std::string ExtractIndexedBases(std::string_view sequence, const std::vector<int>& positions) {
    std::string result;
    result.reserve(positions.size());
    for (const int position : positions) {
        if (position >= 0 && static_cast<std::size_t>(position) < sequence.size()) {
            result.push_back(sequence[static_cast<std::size_t>(position)]);
        }
    }
    return result;
}

inline bool MatchesIupacWithMismatches(
    std::string_view sequence,
    const std::size_t start_offset,
    std::string_view primer,
    const int max_mismatches) {
    if (start_offset + primer.size() > sequence.size()) {
        return false;
    }

    int mismatch_count = 0;
    for (std::size_t index = 0; index < primer.size(); ++index) {
        const char read_base = ToUpper(sequence[start_offset + index]);
        const char primer_base = ToUpper(primer[index]);
        if ((BaseMask(read_base) & BaseMask(primer_base)) == 0) {
            ++mismatch_count;
            if (mismatch_count > max_mismatches) {
                return false;
            }
        }
    }
    return true;
}

inline int FindPrimerFromLeft(std::string_view sequence, std::string_view primer, int max_mismatches) {
    if (primer.size() > sequence.size()) {
        return -1;
    }

    const std::size_t last_start = sequence.size() - primer.size();
    for (std::size_t start = 0; start <= last_start; ++start) {
        if (MatchesIupacWithMismatches(sequence, start, primer, max_mismatches)) {
            return static_cast<int>(start);
        }
    }
    return -1;
}

inline int FindPrimerFromRight(std::string_view sequence, std::string_view primer, int max_mismatches) {
    if (primer.size() > sequence.size()) {
        return -1;
    }

    for (std::size_t start = sequence.size() - primer.size() + 1; start-- > 0;) {
        if (MatchesIupacWithMismatches(sequence, start, primer, max_mismatches)) {
            return static_cast<int>(start);
        }
    }
    return -1;
}

}  // namespace ont::seq
