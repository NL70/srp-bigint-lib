#include "BigIntUtils.hpp"
#include <random>
#include <sstream>
#include <iomanip>
#include <algorithm>

namespace BigIntUtils {

std::string randomHex(size_t bits) {
    size_t hexDigits = (bits + 3) / 4; // 4 bits per hex digit
    std::string s;

    std::mt19937_64 rng(std::random_device{}());
    std::uniform_int_distribution<uint64_t> dist(0, 0xFFFFFFFFFFFFFFFFULL);

    size_t chunks = (hexDigits + 15) / 16;

    for (size_t i = 0; i < chunks; i++) {
        uint64_t val = dist(rng);
        std::stringstream ss;
        ss << std::hex << std::setw(16) << std::setfill('0') << val;
        s = ss.str() + s; // prepend to maintain MSB first
    }

    // Trim extra digits to match exact requested size
    if (s.size() > hexDigits) s = s.substr(s.size() - hexDigits);

    // Make sure first hex digit is non-zero
    if (s[0] == '0') s[0] = '1';

    return s;
}

std::string padHex(const std::string& hexStr, size_t bits) {
    size_t hexDigits = (bits + 3) / 4;
    if (hexStr.size() >= hexDigits) return hexStr;
    return std::string(hexDigits - hexStr.size(), '0') + hexStr;
}

} // namespace BigIntUtils
