#ifndef BIGINT_HPP
#define BIGINT_HPP

#include <vector>
#include <cstdint>
#include <string>

class BigInt
{
public:
    BigInt(); // zero
    BigInt(uint64_t value);
    BigInt(__int128_t value);
    BigInt(int value);
    BigInt(const std::string &hexStr); // arbitrary-length hex

    // Arithmetic
    BigInt operator+(const BigInt &other) const;
    BigInt operator-(const BigInt &other) const;
    BigInt operator*(const BigInt &other) const;
    BigInt operator/(const BigInt &other) const; // long division
    BigInt operator%(const BigInt &other) const;

    // Arithmetic assignment operators
    BigInt &operator+=(const BigInt &other);
    BigInt &operator-=(const BigInt &other);
    BigInt &operator*=(const BigInt &other);
    BigInt &operator/=(const BigInt &other);
    BigInt &operator%=(const BigInt &other);

    BigInt schoolbookMultiply(const BigInt &other) const;

    // Comparison
    bool operator<(const BigInt &other) const;
    bool operator>(const BigInt &other) const;
    bool operator==(const BigInt &other) const;
    bool operator!=(const BigInt &other) const;
    bool operator<=(const BigInt &other) const;
    bool operator>=(const BigInt &other) const;

    // Bit shifting
    BigInt operator<<(size_t bits) const;
    BigInt operator>>(size_t bits) const;

    // Bit shifting assignment operators
    BigInt &operator<<=(size_t bits);
    BigInt &operator>>=(size_t bits);

    // NT funcs
    static BigInt gcd(const BigInt &a, const BigInt &b);
    static BigInt binGCD(const BigInt &a, const BigInt &b);
    static BigInt invert(const BigInt &a, const BigInt &m);

    // Utilities
    std::string toHex() const;
    bool isEven() const;
    size_t bitLength() const;

private:
    std::vector<uint64_t> limbs; // LSB first
    void normalize();

    // Helpers
    static uint64_t hexChunkToUint64(const std::string &chunk);
    static std::pair<BigInt, BigInt> divmod(const BigInt &dividend, const BigInt &divisor);
};

#endif // BIGINT_HPP
