#include <catch2/catch_test_macros.hpp>
#include <catch2/benchmark/catch_benchmark.hpp>
#include <catch2/catch_session.hpp>
#include "naive.hpp"
#include <string>


TEST_CASE("Functional Tests for Naive and Wavefront Alignments") {
    int x = 4, o = 6, e = 2;

    SECTION("Identical Sequences") {
        std::string a = "GATTACA";
        std::string b = "GATTACA";
        REQUIRE(naive(a, b, x, o, e) == 0);
        REQUIRE(wavefront(a, b, x, o, e) == 0);
    }

    SECTION("Completely Different Sequences") {
        std::string a = "GATTACA";
        std::string b = "CTGACGT";
        REQUIRE(naive(a, b, x, o, e) > 0);
        REQUIRE(wavefront(a, b, x, o, e) > 0);
    }

    SECTION("Single Insertion") {
        std::string a = "GATTACA";
        std::string b = "GATTTACA";
        int expected_cost = o;
        REQUIRE(naive(a, b, x, o, e) == expected_cost);
        REQUIRE(wavefront(a, b, x, o, e) == expected_cost);
    }

    SECTION("Single Deletion") {
        std::string a = "GATTACA";
        std::string b = "GATA";
        int expected_cost = o + 2 * e; // Corrected expected cost
        REQUIRE(naive(a, b, x, o, e) == expected_cost);
        REQUIRE(wavefront(a, b, x, o, e) == expected_cost);
    }

    SECTION("Single Mismatch") {
        std::string a = "GATTACA";
        std::string b = "GACTACA";
        int expected_cost = x;
        REQUIRE(naive(a, b, x, o, e) == expected_cost);
        REQUIRE(wavefront(a, b, x, o, e) == expected_cost);
    }

    SECTION("Long Sequence with Single Mismatch") {
        std::string long_a(1000, 'A');
        std::string long_b = long_a;
        long_b[500] = 'T';
        int expected_cost = x;
        REQUIRE(naive(long_a, long_b, x, o, e) == expected_cost);
        REQUIRE(wavefront(long_a, long_b, x, o, e) == expected_cost);
    }

    SECTION("Empty Sequences") {
        std::string a = "";
        std::string b = "";
        REQUIRE(naive(a, b, x, o, e) == 0);
        REQUIRE(wavefront(a, b, x, o, e) == 0);
    }

    SECTION("One Empty Sequence") {
        std::string a = "GATTACA";
        std::string b = "";
        int expected_cost = o + e * (a.size() - 1);
        REQUIRE(naive(a, b, x, o, e) == expected_cost);
        REQUIRE(wavefront(a, b, x, o, e) == expected_cost);
    }
}


// Benchmarking Section for Performance Measurement
TEST_CASE("Benchmarking Naive and Wavefront Alignment Implementations", "[benchmark]") {
    int x = 4, o = 6, e = 2;

    std::string short_a = "GATTACA";
    std::string short_b = "CTGACGT";

    std::string long_a(1000, 'A');
    std::string long_b = long_a;
    long_b[500] = 'T'; // Single mismatch in the middle

    BENCHMARK("Naive Alignment - Short Sequences") {
        return naive(short_a, short_b, x, o, e);
    };

    BENCHMARK("Wavefront Alignment - Short Sequences") {
        return wavefront(short_a, short_b, x, o, e);
    };

    BENCHMARK("Naive Alignment - Long Sequences") {
        return naive(long_a, long_b, x, o, e);
    };

    BENCHMARK("Wavefront Alignment - Long Sequences") {
        return wavefront(long_a, long_b, x, o, e);
    };
}

int main(int argc, char* argv[]) {
    Catch::Session session; // Create a Catch2 session to run the tests

    // Use the session to run all test cases
    int result = session.run(argc, argv);

    return (result < 0xff ? result : 0xff); // Ensure the return value is within range
}