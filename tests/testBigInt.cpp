#include <benchmark/benchmark.h>
#include <gmpxx.h>
#include "BigInt.hpp"
#include "BigIntUtils.hpp"
#include <iostream>

static void BM_BigInt_Add(benchmark::State& state) {
    int BIT_SIZE = state.range(0);
    std::string hexA = BigIntUtils::randomHex(BIT_SIZE);
    std::string hexB = BigIntUtils::randomHex(BIT_SIZE);

    BigInt a(hexA), b(hexB);

    for (auto _ : state) {
        BigInt sum = a + b;
        benchmark::DoNotOptimize(sum);
    }
}

static void BM_GMP_Add(benchmark::State& state) {
    int BIT_SIZE = state.range(0);
    std::string hexA = BigIntUtils::randomHex(BIT_SIZE);
    std::string hexB = BigIntUtils::randomHex(BIT_SIZE);

    mpz_class a("0x" + hexA), b("0x" + hexB);

    for (auto _ : state) {
        mpz_class sum = a + b;
        benchmark::DoNotOptimize(sum);
    }
}

static void BM_BigInt_Minus(benchmark::State& state) {
    int BIT_SIZE = state.range(0);
    std::string hexA = BigIntUtils::randomHex(BIT_SIZE);
    std::string hexB = BigIntUtils::randomHex(BIT_SIZE);

    BigInt a(hexA), b(hexB);

    if (a < b) std::swap(a, b);

    for (auto _ : state) {
        BigInt diff = a - b;
        benchmark::DoNotOptimize(diff);
    }
}

static void BM_GMP_Minus(benchmark::State& state) {
    int BIT_SIZE = state.range(0);
    std::string hexA = BigIntUtils::randomHex(BIT_SIZE);
    std::string hexB = BigIntUtils::randomHex(BIT_SIZE);

    mpz_class a("0x" + hexA), b("0x" + hexB);

    if (a < b) std::swap(a, b);

    for (auto _ : state) {
        mpz_class diff = a - b;
        benchmark::DoNotOptimize(diff);
    }
}

static void BM_BigInt_Karatsuba_Multiply(benchmark::State& state) {
    int BIT_SIZE = state.range(0);
    std::string hexA = BigIntUtils::randomHex(BIT_SIZE);
    std::string hexB = BigIntUtils::randomHex(BIT_SIZE);

    BigInt a(hexA), b(hexB);

    for (auto _ : state) {
        BigInt prod = a * b;
        benchmark::DoNotOptimize(prod);
    }
}

static void BM_BigInt_Schoolbook_Multiply(benchmark::State& state) {
    int BIT_SIZE = state.range(0);
    std::string hexA = BigIntUtils::randomHex(BIT_SIZE);
    std::string hexB = BigIntUtils::randomHex(BIT_SIZE);

    BigInt a(hexA), b(hexB);

    for (auto _ : state) {
        BigInt prod = a.schoolbookMultiply(b);
        benchmark::DoNotOptimize(prod);
    }
}

static void BM_GMP_Multiply(benchmark::State& state) {
    int BIT_SIZE = state.range(0);
    std::string hexA = BigIntUtils::randomHex(BIT_SIZE);
    std::string hexB = BigIntUtils::randomHex(BIT_SIZE);

    mpz_class a("0x" + hexA), b("0x" + hexB);

    for (auto _ : state) {
        mpz_class prod = a * b;
        benchmark::DoNotOptimize(prod);
    }
}

static void BM_BigInt_Right_Shift(benchmark::State& state) {
    int BIT_SIZE = state.range(0);
    std::string hexA = BigIntUtils::randomHex(BIT_SIZE);
    BigInt a(hexA);

    size_t shiftAmount = (std::rand() % 512) + 1;

    for (auto _ : state) {
        BigInt rightShift = a >> shiftAmount;
        benchmark::DoNotOptimize(rightShift);
    }
}

static void BM_GMP_Right_Shift(benchmark::State& state) {
    int BIT_SIZE = state.range(0);
    std::string hexA = BigIntUtils::randomHex(BIT_SIZE);
    mpz_class a("0x" + hexA);

    size_t shiftAmount = (std::rand() % 512) + 1;

    for (auto _ : state) {
        mpz_class rightShift = a >> shiftAmount;
        benchmark::DoNotOptimize(rightShift);
    }
}

static void BM_BigInt_Left_Shift(benchmark::State& state) {
    int BIT_SIZE = state.range(0);
    std::string hexA = BigIntUtils::randomHex(BIT_SIZE);
    BigInt a(hexA);

    size_t shiftAmount = (std::rand() % 512) + 1;

    for (auto _ : state) {
        BigInt leftShift = a << shiftAmount;
        benchmark::DoNotOptimize(leftShift);
    }
}

static void BM_GMP_Left_Shift(benchmark::State& state) {
    int BIT_SIZE = state.range(0);
    std::string hexA = BigIntUtils::randomHex(BIT_SIZE);
    mpz_class a("0x" + hexA);

    size_t shiftAmount = (std::rand() % 512) + 1;

    for (auto _ : state) {
        mpz_class leftShift = a << shiftAmount;
        benchmark::DoNotOptimize(leftShift);
    }
}

static void BM_BigInt_Divide(benchmark::State& state) {
    int BIT_SIZE = state.range(0);
    std::string hexA = BigIntUtils::randomHex(BIT_SIZE * 2);
    std::string hexB = BigIntUtils::randomHex(BIT_SIZE);

    BigInt a(hexA), b(hexB);

    for (auto _ : state) {
        BigInt quot = a / b;
        benchmark::DoNotOptimize(quot);
    }
}

static void BM_GMP_Divide(benchmark::State& state) {
    int BIT_SIZE = state.range(0);
    std::string hexA = BigIntUtils::randomHex(BIT_SIZE * 2);
    std::string hexB = BigIntUtils::randomHex(BIT_SIZE);

    mpz_class a("0x" + hexA), b("0x" + hexB);

    for (auto _ : state) {
        mpz_class quot = a / b;
        benchmark::DoNotOptimize(quot);
    }
}

static void BM_BigInt_Reduction(benchmark::State& state) {
    int BIT_SIZE = state.range(0);
    std::string hexA = BigIntUtils::randomHex(BIT_SIZE * 2);
    std::string hexB = BigIntUtils::randomHex(BIT_SIZE);

    BigInt a(hexA), b(hexB);

    for (auto _ : state) {
        BigInt mod = a % b;
        benchmark::DoNotOptimize(mod);
    }
}

static void BM_GMP_Reduction(benchmark::State& state) {
    int BIT_SIZE = state.range(0);
    std::string hexA = BigIntUtils::randomHex(BIT_SIZE * 2);
    std::string hexB = BigIntUtils::randomHex(BIT_SIZE);

    mpz_class a("0x" + hexA), b("0x" + hexB);

    for (auto _ : state) {
        mpz_class mod = a % b;
        benchmark::DoNotOptimize(mod);
    }
}

static void BM_BigInt_Lehmer_GCD(benchmark::State& state) {
    int BIT_SIZE = state.range(0);
    std::string hexA = BigIntUtils::randomHex(BIT_SIZE);
    std::string hexB = BigIntUtils::randomHex(BIT_SIZE);

    BigInt a(hexA), b(hexB);

    for (auto _ : state) {
        BigInt mod = BigInt::gcd(a, b);
        benchmark::DoNotOptimize(mod);
    }
}

static void BM_BigInt_Binary_GCD(benchmark::State& state) {
    int BIT_SIZE = state.range(0);
    std::string hexA = BigIntUtils::randomHex(BIT_SIZE);
    std::string hexB = BigIntUtils::randomHex(BIT_SIZE);

    BigInt a(hexA), b(hexB);

    for (auto _ : state) {
        BigInt GCD = BigInt::binGCD(a, b);
        benchmark::DoNotOptimize(GCD);
    }
}

static void BM_GMP_GCD(benchmark::State& state) {
    int BIT_SIZE = state.range(0);
    std::string hexA = BigIntUtils::randomHex(BIT_SIZE );
    std::string hexB = BigIntUtils::randomHex(BIT_SIZE);

    mpz_class a("0x" + hexA), b("0x" + hexB);

    for (auto _ : state) {
        mpz_class GCD = gcd(a, b);
        benchmark::DoNotOptimize(GCD);
    }
}

static void BM_BigInt_Invert(benchmark::State& state) {
    int BIT_SIZE = state.range(0);
    BigInt a, b;

    while (true) {
        std::string hexA = BigIntUtils::randomHex(BIT_SIZE);
        std::string hexB = BigIntUtils::randomHex(BIT_SIZE);

        a = BigInt(hexA);
        b = BigInt(hexB);
        if (BigInt::gcd(a, b) == 1) break;
    }

    for (auto _ : state) {
        BigInt inverse = BigInt::invert(a, b);
        benchmark::DoNotOptimize(inverse);
    }
}

static void BM_GMP_Invert(benchmark::State& state) {
    int BIT_SIZE = state.range(0);
    mpz_class a, b;

    while (true) {
        std::string hexA = BigIntUtils::randomHex(BIT_SIZE);
        std::string hexB = BigIntUtils::randomHex(BIT_SIZE);

        a = mpz_class("0x" + hexA);
        b = mpz_class("0x" + hexB);
        if (gcd(a, b) == 1) break;
    }

    for (auto _ : state) {
        mpz_class gmpInvert;
        int success = mpz_invert(gmpInvert.get_mpz_t(), a.get_mpz_t(), b.get_mpz_t());
        benchmark::DoNotOptimize(gmpInvert);
        benchmark::DoNotOptimize(success);
    }
}

BENCHMARK(BM_BigInt_Karatsuba_Multiply)
    ->DenseRange(1<<10, 1<<14, 1<<10);   

BENCHMARK(BM_GMP_Multiply)
    ->DenseRange(1<<10, 1<<14, 1<<10);   

BENCHMARK_MAIN();
