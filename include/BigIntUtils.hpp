#ifndef BIGINT_UTILS_HPP
#define BIGINT_UTILS_HPP

#include <string>
#include <cstddef>

namespace BigIntUtils {

    // Generate a random hex string representing a number with `bits` bits
    std::string randomHex(size_t bits);

    // Optionally, pad a hex string to a specific bit length
    std::string padHex(const std::string& hexStr, size_t bits);

}

#endif // BIGINT_UTILS_HPP
