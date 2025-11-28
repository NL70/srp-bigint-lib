#include "BigInt.hpp"
#include <sstream>
#include <iomanip>
#include <algorithm>
#include <stdexcept>
#include <cctype>
#include <iostream>

// ---------------- Constructors ----------------
BigInt::BigInt() { limbs.push_back(0); }

BigInt::BigInt(uint64_t value) { limbs.push_back(value); }

BigInt::BigInt(__int128_t value)
{
    // Implementation to create a BigInt from a 128-bit integer
    // It might have one or two limbs

    if (value < 0)
    {
        throw std::invalid_argument("Negative value cannot be assigned to unsigned BigInt.");
    }
    limbs.clear();
    if (value == 0)
        return;
    uint64_t lower = static_cast<uint64_t>(value);
    uint64_t upper = static_cast<uint64_t>(value >> 64);
    limbs.push_back(lower);
    if (upper > 0)
    {
        limbs.push_back(upper);
    }
}

BigInt::BigInt(int value) : BigInt(static_cast<uint64_t>(value)) {}

BigInt::BigInt(const std::string &hexStr)
{
    std::string s = hexStr;
    if (s.rfind("0x", 0) == 0 || s.rfind("0X", 0) == 0)
        s = s.substr(2);

    size_t len = s.size();
    if (len == 0)
    {
        limbs.push_back(0);
        return;
    }
    size_t rem = len % 16;
    if (rem != 0)
        s = std::string(16 - rem, '0') + s;

    for (ssize_t pos = (ssize_t)s.size() - 16; pos >= 0; pos -= 16)
    {
        std::string chunk = s.substr(pos, 16);
        limbs.push_back(hexChunkToUint64(chunk));
        if (pos == 0)
            break;
    }

    normalize();
}

// ---------------- Helper ----------------
void BigInt::normalize()
{
    while (limbs.size() > 1 && limbs.back() == 0)
        limbs.pop_back();
}

uint64_t BigInt::hexChunkToUint64(const std::string &chunk)
{
    return std::stoull(chunk, nullptr, 16);
}

// ---------------- Arithmetic ----------------
BigInt &BigInt::operator+=(const BigInt &other)
{
    size_t n = std::max(limbs.size(), other.limbs.size());
    limbs.resize(n+1, 0); // preallocate an extra limb 
    uint64_t carry = 0;
    for (size_t i = 0; i < n; i++)
    {
        uint64_t b = i < other.limbs.size() ? other.limbs[i] : 0;
        __uint128_t sum = (__uint128_t)limbs[i] + b + carry;
        limbs[i] = (uint64_t)sum;
        carry = (uint64_t)(sum >> 64);
    }
    if (carry) {
        limbs[n] = carry;
    } else
        normalize();
    return *this;
}

BigInt &BigInt::operator-=(const BigInt &other)
{
    if (*this < other)
        throw std::runtime_error("Negative result not supported");

    uint64_t borrow = 0;
    for (size_t i = 0; i < limbs.size(); i++)
    {
        uint64_t a = limbs[i];
        uint64_t b = i < other.limbs.size() ? other.limbs[i] : 0;
        __uint128_t sub = (__uint128_t)a - b - borrow;
        limbs[i] = (uint64_t)sub;
        borrow = (sub >> 127) & 1;
    }
    normalize();
    return *this;
}

BigInt &BigInt::operator*=(const BigInt &other)
{
    *this = *this * other;
    return *this;
}

BigInt &BigInt::operator/=(const BigInt &other)
{
    *this = divmod(*this, other).first;
    return *this;
}

BigInt &BigInt::operator%=(const BigInt &other)
{
    *this = divmod(*this, other).second;
    return *this;
}

BigInt BigInt::operator+(const BigInt &other) const
{
    BigInt result = *this;
    result += other;
    return result;
}

BigInt BigInt::operator-(const BigInt &other) const
{
    BigInt result = *this;
    result -= other;
    return result;
}

BigInt BigInt::operator*(const BigInt &other) const
{
    if (this->limbs.size() <= 64 || other.limbs.size() <= 64)
    {
        return schoolbookMultiply(other);
    }

    size_t n = std::max(this->limbs.size(), other.limbs.size());
    size_t m = (n + 1) / 2;

    BigInt low1, high1;
    if (limbs.size() > m)
    {
        low1.limbs.assign(limbs.begin(), limbs.begin() + m);
        high1.limbs.assign(limbs.begin() + m, limbs.end());
    }
    else
    {
        low1 = *this;
        high1 = BigInt(0);
    }

    BigInt low2, high2;
    if (other.limbs.size() > m)
    {
        low2.limbs.assign(other.limbs.begin(), other.limbs.begin() + m);
        high2.limbs.assign(other.limbs.begin() + m, other.limbs.end());
    }
    else
    {
        low2 = other;
        high2 = BigInt(0);
    }

    BigInt z0 = low1 * low2;
    BigInt z2 = high1 * high2;

    BigInt s1 = low1;
    s1 += high1;
    BigInt s2 = low2;
    s2 += high2;

    BigInt z1 = s1 * s2;
    z1 -= z0;
    z1 -= z2;

    z2 <<= (2 * m * 64);
    z1 <<= (m * 64);

    BigInt result = z0;
    result += z1;
    result += z2;

    return result;
}

BigInt BigInt::operator/(const BigInt &divisor) const
{
    return divmod(*this, divisor).first;
}

BigInt BigInt::operator%(const BigInt &divisor) const
{
    return divmod(*this, divisor).second;
}

// ---------------- Comparison ----------------
bool BigInt::operator<(const BigInt &other) const
{
    if (limbs.size() != other.limbs.size())
        return limbs.size() < other.limbs.size();
    for (ssize_t i = limbs.size() - 1; i >= 0; i--)
    {
        if (limbs[i] != other.limbs[i])
            return limbs[i] < other.limbs[i];
    }
    return false;
}

bool BigInt::operator>(const BigInt &other) const
{
    if (limbs.size() != other.limbs.size())
        return limbs.size() > other.limbs.size();
    for (ssize_t i = limbs.size() - 1; i >= 0; i--)
    {
        if (limbs[i] != other.limbs[i])
            return limbs[i] > other.limbs[i];
    }
    return false;
}

bool BigInt::operator==(const BigInt &other) const { return limbs == other.limbs; }
bool BigInt::operator!=(const BigInt &other) const { return !(*this == other); }
bool BigInt::operator<=(const BigInt &other) const { return !(*this > other); }
bool BigInt::operator>=(const BigInt &other) const { return !(*this < other); }

// ---------------- Bit shifts ----------------
BigInt &BigInt::operator<<=(size_t bits)
{
    if (limbs.size() == 1 && limbs[0] == 0 || bits == 0)
    {
        return *this;
    }

    size_t limbShift = bits / 64;
    size_t bitShift = bits % 64;

    std::vector<uint64_t> result_limbs;
    result_limbs.resize(limbs.size() + limbShift + 1, 0);

    uint64_t carry = 0;
    for (size_t i = 0; i < limbs.size(); ++i)
    {
        __uint128_t shifted = (__uint128_t)limbs[i] << bitShift;
        result_limbs[i + limbShift] = (uint64_t)(shifted & 0xFFFFFFFFFFFFFFFFULL) | carry;
        carry = (uint64_t)(shifted >> 64);
    }

    if (carry)
        result_limbs[limbs.size() + limbShift] = carry;

    limbs = std::move(result_limbs);
    normalize();
    return *this;
}

BigInt &BigInt::operator>>=(size_t bits)
{
    if (limbs.size() == 1 && limbs[0] == 0 || bits == 0)
    {
        return *this;
    }

    size_t limbShift = bits / 64;
    size_t bitShift = bits % 64;

    if (limbShift >= limbs.size())
    {
        limbs = {0};
        return *this;
    }

    std::vector<uint64_t> result_limbs;
    result_limbs.resize(limbs.size() - limbShift);

    uint64_t carry = 0;
    for (ssize_t i = (ssize_t)limbs.size() - 1; i >= (ssize_t)limbShift; --i)
    {
        uint64_t current = limbs[i];
        result_limbs[i - limbShift] = (current >> bitShift) | (carry << (64 - bitShift));
        carry = current & ((1ULL << bitShift) - 1ULL);
    }

    limbs = std::move(result_limbs);
    normalize();
    return *this;
}

BigInt BigInt::operator<<(size_t bits) const
{
    BigInt result = *this;
    result <<= bits;
    return result;
}

BigInt BigInt::operator>>(size_t bits) const
{
    BigInt result = *this;
    result >>= bits;
    return result;
}

// ---------------- NT funcs ----------------

BigInt BigInt::gcd(const BigInt &a, const BigInt &b)
{
    BigInt u = a, v = b;

    if (u < v)
        std::swap(u, v);

    if (v == BigInt(0))
        return u;

    BigInt r, term1, term2, next_u, next_v;
    BigInt A_big, B_big, C_big, D_big;

    while (true)
    {

        if (v.limbs.size() <= 2 || u.limbs.size() != v.limbs.size())
        {
            r = u % v;
            u = v;
            v = r;
            if (v == BigInt(0))
                return u;
            continue;
        }

        // Combine top two limbs into a uint128
        unsigned __int128 u_top = (static_cast<unsigned __int128>(u.limbs.back()) << 64) | u.limbs[u.limbs.size() - 2];
        unsigned __int128 v_top = (static_cast<unsigned __int128>(v.limbs.back()) << 64) | v.limbs[v.limbs.size() - 2];

        // Get top 127 bits so it fits in int128
        __int128_t x = static_cast<__int128_t>(u_top >> 1);
        __int128_t y = static_cast<__int128_t>(v_top >> 1);

        if (y == 0)
        {
            r = u % v;
            u = v;
            v = r;
            if (v == BigInt(0))
                return u;
            continue;
        }

        __int128_t A = 1, B = 0;
        __int128_t C = 0, D = 1;

        while (true)
        {
            if ((y + C) == 0 || (y + D) == 0)
                break;

            __int128_t w1 = (x + A) / (y + C);
            __int128_t w2 = (x + B) / (y + D);

            if (w1 != w2)
                break;

            __int128_t w = w1;
            __int128_t temp_A = A;
            A = C;
            C = temp_A - w * C;
            __int128_t temp_B = B;
            B = D;
            D = temp_B - w * D;
            __int128_t temp_x = x;
            x = y;
            y = temp_x - w * y;
        }

        // if B is 0, the inner loop made no progress
        if (B == 0)
        {
            r = u % v;
            u = v;
            v = r;
            if (v == BigInt(0))
                return u;
            continue;
        }

        // apply (u, v) * M
        A_big = BigInt((A < 0) ? -A : A);
        B_big = BigInt((B < 0) ? -B : B);
        C_big = BigInt((C < 0) ? -C : C);
        D_big = BigInt((D < 0) ? -D : D);

        term1 = u * A_big;
        term2 = v * B_big;
        if ((A >= 0) == (B >= 0))
            next_u = term1 + term2;
        else if (term1 > term2)
            next_u = term1 - term2;
        else
            next_u = term2 - term1;

        term1 = u * C_big;
        term2 = v * D_big;
        if ((C >= 0) == (D >= 0))
            next_v = term1 + term2;
        else if (term1 > term2)
            next_v = term1 - term2;
        else
            next_v = term2 - term1;

        u = next_u;
        v = next_v;

        if (v == BigInt(0))
        {
            return u;
        }
    }
}

BigInt BigInt::binGCD(const BigInt &a, const BigInt &b)
{
    BigInt u = a;
    BigInt v = b;

    if (u == BigInt(0))
        return v;

    if (v == BigInt(0))
        return u;

    int k = 0;

    while (u.isEven() && v.isEven())
    {
        u >>= 1;
        v >>= 1;
        k++;
    }

    while (u.isEven())
        u >>= 1;

    while (v != BigInt(0))
    {
        while (v.isEven())
            v >>= 1;
        if (u > v)
            std::swap(u, v);
        v -= u;
    }

    u <<= k;

    return u;
}

//  (a - b) mod m
BigInt mod_sub(const BigInt &a, const BigInt &b, const BigInt &m)
{
    if (a >= b)
    {
        return a - b;
    }
    else
    {
        return m - (b - a);
    }
}

BigInt BigInt::invert(const BigInt &a, const BigInt &m)
{
    // ax + my = gcd(a, m) = 1
    BigInt u = a, v = m;
    BigInt x_u = BigInt(1), x_v = BigInt(0);

    BigInt q, r, temp_x;
    BigInt A_big, B_big, C_big, D_big;
    BigInt term1, term2, next_u, next_v, next_x_u, next_x_v;

    while (v != BigInt(0))
    {
        if (v.limbs.size() <= 2 || u.limbs.size() != v.limbs.size())
        {
            auto [q, r] = BigInt::divmod(u, v);
            u = v;
            v = r;

            // x_new = x_old - q * x_prev (mod m)
            term1 = (q * x_v) % m;
            temp_x = mod_sub(x_u, term1, m);
            x_u = x_v;
            x_v = temp_x;

            continue;
        }

        unsigned __int128 u_top = (static_cast<unsigned __int128>(u.limbs.back()) << 64) | u.limbs[u.limbs.size() - 2];
        unsigned __int128 v_top = (static_cast<unsigned __int128>(v.limbs.back()) << 64) | v.limbs[v.limbs.size() - 2];

        __int128_t x = static_cast<__int128_t>(u_top >> 1);
        __int128_t y = static_cast<__int128_t>(v_top >> 1);

        if (y == 0)
        {
            auto [q, r] = BigInt::divmod(u, v);
            u = v;
            v = r;
            term1 = (q * x_v) % m;
            temp_x = mod_sub(x_u, term1, m);
            x_u = x_v;
            x_v = temp_x;
            continue;
        }

        __int128_t A = 1, B = 0;
        __int128_t C = 0, D = 1;

        while (true)
        {
            if ((y + C) == 0 || (y + D) == 0)
                break;
            __int128_t w1 = (x + A) / (y + C), w2 = (x + B) / (y + D);
            if (w1 != w2)
                break;
            __int128_t w = w1;
            __int128_t temp_A = A;
            A = C;
            C = temp_A - w * C;
            __int128_t temp_B = B;
            B = D;
            D = temp_B - w * D;
            __int128_t temp_x = x;
            x = y;
            y = temp_x - w * y;
        }

        if (B == 0)
        {
            auto [q, r] = BigInt::divmod(u, v);
            u = v;
            v = r;
            term1 = (q * x_v) % m;
            temp_x = mod_sub(x_u, term1, m);
            x_u = x_v;
            x_v = temp_x;
            continue;
        }

        A_big = BigInt((A < 0) ? -A : A);
        B_big = BigInt((B < 0) ? -B : B);
        C_big = BigInt((C < 0) ? -C : C);
        D_big = BigInt((D < 0) ? -D : D);

        // Calculate  x_new = (x_u*A + x_v*B) mod m
        term1 = (x_u * A_big) % m;
        term2 = (x_v * B_big) % m;
        if ((A >= 0) == (B >= 0))
        {
            next_x_u = (term1 + term2) % m;
            if (next_x_u >= m)
            {
                next_x_u -= m;
            }
        }
        else if (A >= 0)
            next_x_u = mod_sub(term1, term2, m); // A is pos, B is neg
        else
            next_x_u = mod_sub(term2, term1, m); // A is neg, B is pos

        term1 = (x_u * C_big) % m;
        term2 = (x_v * D_big) % m;
        if ((C >= 0) == (D >= 0))
        {
            next_x_u = (term1 + term2) % m;
            if (next_x_u >= m)
            {
                next_x_u -= m;
            }
        }
        else if (C >= 0)
            next_x_v = mod_sub(term1, term2, m); // C is pos, D is neg
        else
            next_x_v = mod_sub(term2, term1, m); // C is neg, D is pos

        term1 = u * A_big;
        term2 = v * B_big;
        if ((A >= 0) == (B >= 0))
            next_u = term1 + term2;
        else if (term1 > term2)
            next_u = term1 - term2;
        else
            next_u = term2 - term1;

        term1 = u * C_big;
        term2 = v * D_big;
        if ((C >= 0) == (D >= 0))
            next_v = term1 + term2;
        else if (term1 > term2)
            next_v = term1 - term2;
        else
            next_v = term2 - term1;

        u = next_u;
        v = next_v;
        x_u = next_x_u;
        x_v = next_x_v;
    }

    if (u != BigInt(1))
    {
        return BigInt(0);
    }

    return x_u;
}

// ---------------- Utilities ----------------
std::string BigInt::toHex() const
{
    if (limbs.empty() || (limbs.size() == 1 && limbs[0] == 0))
        return "0";

    std::stringstream ss;
    ss << std::hex;

    ss << limbs.back();

    for (size_t i = limbs.size() - 1; i-- > 0;)
    {
        ss << std::setw(16) << std::setfill('0') << limbs[i];
    }

    return ss.str();
}

bool BigInt::isEven() const
{
    if (limbs.empty() || (limbs.size() == 1 && limbs[0] == 0))
        return true;
    return (limbs[0] & 1) == 0;
}

size_t BigInt::bitLength() const
{
    if (limbs.empty() || (limbs.size() == 1 && limbs[0] == 0))
        return 0;
    uint64_t top = limbs.back();
    size_t bits = 0;
    while (top)
    {
        top >>= 1;
        bits++;
    }
    return bits + 64 * (limbs.size() - 1);
}

// ---------------- Helpers ----------------
BigInt BigInt::schoolbookMultiply(const BigInt &other) const
{
    BigInt result;
    result.limbs.assign(limbs.size() + other.limbs.size(), 0);

    for (size_t i = 0; i < limbs.size(); i++)
    {
        uint64_t carry = 0;
        for (size_t j = 0; j < other.limbs.size(); j++)
        {
            __uint128_t prod = (__uint128_t)limbs[i] * other.limbs[j] + result.limbs[i + j] + carry;
            result.limbs[i + j] = (uint64_t)prod;
            carry = (uint64_t)(prod >> 64);
        }
        result.limbs[i + other.limbs.size()] = carry;
    }

    result.normalize();
    return result;
}
std::pair<BigInt, BigInt> BigInt::divmod(const BigInt &dividend, const BigInt &divisor)
{
    if (divisor == BigInt(0))
        throw std::runtime_error("Division by zero");
    if (dividend < divisor)
        return {BigInt(0), dividend};
    if (divisor.limbs.size() == 1 && divisor.limbs[0] == 1)
        return {dividend, BigInt(0)};

    BigInt u = dividend;
    BigInt v = divisor;

    int n = v.limbs.size();
    int m = u.limbs.size() - n;

    int shift = 0;
    if (v.limbs.back() < (1ULL << 63))
    { // Check if MSB is not set
        shift = __builtin_clzll(v.limbs.back());
        u <<= shift;
        v <<= shift;
    }

    if (u.limbs.size() == (size_t)m + n)
        u.limbs.push_back(0);

    BigInt Q;
    Q.limbs.assign(m + 1, 0);

    const unsigned __int128 BASE = (unsigned __int128)1 << 64;

    for (int i = m; i >= 0; i--)
    {
        unsigned __int128 numer = ((unsigned __int128)u.limbs[i + n] << 64) | u.limbs[i + n - 1];
        unsigned __int128 qhat;

        if (v.limbs.size() > (size_t)n - 1)
        { // Guard against access
            qhat = numer / v.limbs[n - 1];
        }
        else
        {
            qhat = 0;
        }

        if (qhat >= BASE)
            qhat = BASE - 1;

        // multiply and subtract
        BigInt temp_sub;
        temp_sub.limbs.assign(n + 1, 0);
        uint64_t k = 0;
        for (int j = 0; j < n; ++j)
        {
            unsigned __int128 p = (unsigned __int128)v.limbs[j] * (uint64_t)qhat + k;
            temp_sub.limbs[j] = p;
            k = p >> 64;
        }
        temp_sub.limbs[n] = k;

        temp_sub.normalize();

        BigInt u_slice;
        u_slice.limbs.assign(u.limbs.begin() + i, u.limbs.begin() + i + n + 1);

        while (u_slice < temp_sub)
        {
            qhat--;
            temp_sub -= v;
        }

        u_slice -= temp_sub;
        Q.limbs[i] = (uint64_t)qhat;

        for (size_t j = 0; j < u_slice.limbs.size(); ++j)
        {
            if (i + j < u.limbs.size())
                u.limbs[i + j] = u_slice.limbs[j];
        }
    }

    BigInt R;
    R.limbs.assign(u.limbs.begin(), u.limbs.begin() + n);
    if (shift > 0)
        R >>= shift;

    Q.normalize();
    R.normalize();
    return {Q, R};
}