// Sequence.hpp
// Header-only nucleotide utilities: ASCII base handling, IUPAC bit-mask
// lookup for ambiguity-aware matching, reverse-complement, and primer
// search helpers. Designed for hot paths in PrimerTrimmer and the DMS
// filter where every read calls these multiple times, so tables are built
// once (Meyers singletons) and reused.
#pragma once

#include <algorithm>
#include <array>
#include <cctype>
#include <cstdint>
#include <string>
#include <string_view>
#include <vector>

namespace ont::seq {

// ASCII-safe uppercase that avoids std::toupper's signed-char UB.
inline char ToUpper(char value) {
    return static_cast<char>(std::toupper(static_cast<unsigned char>(value)));
}

// True for A/C/G/T (case-insensitive); ambiguous IUPAC codes return false.
inline bool IsCanonicalBase(char value) {
    const char base = ToUpper(value);
    return base == 'A' || base == 'C' || base == 'G' || base == 'T';
}

namespace detail {

// 256-entry lookup keyed by ASCII code: each entry is a 4-bit mask over
// {A=1, C=2, G=4, T=8} encoding which canonical bases an IUPAC code allows.
// Built once lazily so repeated calls on the hot path are table lookups.
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

// 256-entry IUPAC complement lookup. Unset entries default to 'N' so any
// stray non-nucleotide byte reverse-complements to an unambiguous no-call
// rather than silently returning 0.
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

// Return the IUPAC base mask for a character (see BaseMaskTable).
inline std::uint8_t BaseMask(char value) {
    return detail::BaseMaskTable()[static_cast<unsigned char>(value)];
}

// Return the IUPAC complement of a base. Unknown inputs map to 'N'.
inline char Complement(char value) {
    return detail::ComplementTable()[static_cast<unsigned char>(value)];
}

// Uppercase a sequence, returning a new string. Used to normalise reads once
// before a primer search instead of upper-casing on every compare.
inline std::string ToUpperCopy(std::string value) {
    std::transform(value.begin(), value.end(), value.begin(), ToUpper);
    return value;
}

// In-place uppercase variant; callers keep ownership of the buffer.
inline void UppercaseInPlace(std::string* value) {
    std::transform(value->begin(), value->end(), value->begin(), ToUpper);
}

// Reverse a string without complementing - used to keep quality strings in
// sync when a read is reverse-complemented.
inline std::string ReverseCopy(std::string_view value) {
    return std::string(value.rbegin(), value.rend());
}

// Alias for ReverseCopy that documents intent at call sites.
inline std::string ReverseQuality(std::string_view quality) {
    return ReverseCopy(quality);
}

// Return the reverse complement of a sequence. IUPAC ambiguous bases are
// complemented via the lookup table (e.g. R<->Y, W<->W).
inline std::string ReverseComplement(std::string_view sequence) {
    std::string reverse_complement;
    reverse_complement.reserve(sequence.size());
    for (auto iterator = sequence.rbegin(); iterator != sequence.rend(); ++iterator) {
        reverse_complement.push_back(Complement(*iterator));
    }
    return reverse_complement;
}

// List the 0-based offsets at which a primer carries an ambiguous IUPAC
// code. Used by UMI extraction to know which bases in the aligned primer
// should be collected as UMI payload vs. match check.
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

// Pick the bases at the given indices out of sequence and concatenate them.
// Out-of-range indices are silently skipped so a short read cannot abort
// extraction mid-way.
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

// Check whether primer aligns to sequence at start_offset with at most
// max_mismatches. IUPAC compatibility uses bit-wise AND of the two masks:
// non-zero means at least one shared base, which is treated as a match.
// Returns false early once the mismatch budget is exceeded.
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

// Scan left-to-right and return the first start offset where primer matches
// sequence within the mismatch budget, or -1 if nothing fits. Used to locate
// the forward primer at the 5' end of a nanopore read.
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

// Scan right-to-left and return the right-most matching start offset, or -1
// if the primer does not fit. Used to anchor the reverse-complement of the
// reverse primer near the 3' end so the trimmed span brackets the insert.
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
