
#include <catch2/catch_test_macros.hpp>
#include <catch2/benchmark/catch_benchmark.hpp>
#include <catch2/catch_session.hpp>

#include "include/naive.hpp"
#include <string>
#include <algorithm>

// Test Suite for Alignment Functions
TEST_CASE("Comprehensive Tests for wfa::naive and wfa::wavefront_dp Alignments") {
    int x = 4, o = 6, e = 2; // Default penalty values

    SECTION("Identical Sequences") {
        std::string a = "GATTACA";
        std::string b = "GATTACA";
        int expected_cost = 0; // 7M
        REQUIRE(wfa::naive(a, b, x, o, e) == expected_cost);
        REQUIRE(wfa::wavefront_dp(a, b, x, o, e) == expected_cost);
    }

    SECTION("Completely Different Sequences") {
        std::string a = "GATTACA";
        std::string b = "CTGACGT";
        int expected_cost = -28; // 7X
        REQUIRE(wfa::naive(a, b, x, o, e) == expected_cost);
        REQUIRE(wfa::wavefront_dp(a, b, x, o, e) == expected_cost);
    }

    SECTION("Single Insertion") {
        std::string a = "GATTACA";
        std::string b = "GATTTACA";
        int expected_cost = -8; // 1I
        REQUIRE(wfa::naive(a, b, x, o, e) == expected_cost);
        REQUIRE(wfa::wavefront_dp(a, b, x, o, e) == expected_cost);
    }

    SECTION("Single Deletion") {
        std::string a = "GATTACA";
        std::string b = "GATA";
        int expected_cost = -12; // 3D
        REQUIRE(wfa::naive(a, b, x, o, e) == expected_cost);
        REQUIRE(wfa::wavefront_dp(a, b, x, o, e) == expected_cost);
    }

    SECTION("Single Mismatch") {
        std::string a = "GATTACA";
        std::string b = "GACTACA";
        int expected_cost = -4; // 1X
        REQUIRE(wfa::naive(a, b, x, o, e) == expected_cost);
        REQUIRE(wfa::wavefront_dp(a, b, x, o, e) == expected_cost);
    }

    SECTION("Long Sequence with Single Mismatch") {
        std::string long_a(1000, 'A');
        std::string long_b = long_a;
        long_b[500] = 'T';
        int expected_cost = -4; // 1X
        REQUIRE(wfa::naive(long_a, long_b, x, o, e) == expected_cost);
        REQUIRE(wfa::wavefront_dp(long_a, long_b, x, o, e) == expected_cost);
    }

    SECTION("Empty Sequences") {
        std::string a = "";
        std::string b = "";
        int expected_cost = 0;
        REQUIRE(wfa::naive(a, b, x, o, e) == expected_cost);
        REQUIRE(wfa::wavefront_dp(a, b, x, o, e) == expected_cost);
    }

    SECTION("One Empty Sequence") {
        std::string a = "GATTACA";
        std::string b = "";
        int expected_cost = -20; // 7D
        REQUIRE(wfa::wavefront_dp(a, b, x, o, e) == expected_cost);
    }

    SECTION("Multiple Gaps and Mismatches") {
        std::string a = "AGCTTAC";
        std::string b = "CGTTAGC";
        int expected_cost = -16; // 1X1X2X
        REQUIRE(wfa::naive(a, b, x, o, e) == expected_cost);
        REQUIRE(wfa::wavefront_dp(a, b, x, o, e) == expected_cost);
    }

    SECTION("Consecutive Gaps") {
        std::string a = "AAAAAA";
        std::string b = "AA";
        int expected_cost = -14; // 4D
        REQUIRE(wfa::naive(a, b, x, o, e) == expected_cost);
        REQUIRE(wfa::wavefront_dp(a, b, x, o, e) == expected_cost);
    }

    SECTION("Alternating Matches and Mismatches") {
        std::string a = "ACACACAC";
        std::string b = "AAAA";
        int expected_cost = -22; // 1X4D1X
        REQUIRE(wfa::naive(a, b, x, o, e) == expected_cost);
        REQUIRE(wfa::wavefront_dp(a, b, x, o, e) == expected_cost);
    }

    SECTION("Long Identical Sequences") {
        std::string a(10000, 'G');
        std::string b(10000, 'G');
        int expected_cost = 0;
        REQUIRE(wfa::naive(a, b, x, o, e) == expected_cost);
        REQUIRE(wfa::wavefront_dp(a, b, x, o, e) == expected_cost);
    }

    SECTION("Repeated Characters with Gaps") {
        std::string a = "TTTTTTTT";
        std::string b = "TTTTGGGG";
        int expected_cost = -16; // 4X
        REQUIRE(wfa::naive(a, b, x, o, e) == expected_cost);
        REQUIRE(wfa::wavefront_dp(a, b, x, o, e) == expected_cost);
    }

    SECTION("Subsequence Alignment") {
        std::string a = "AGCTAGC";
        std::string b = "GCTA";
        int expected_cost = -18; // 1D2D
        REQUIRE(wfa::naive(a, b, x, o, e) == expected_cost);
        REQUIRE(wfa::wavefront_dp(a, b, x, o, e) == expected_cost);
    }

    SECTION("High Mismatch Penalty Preference for Gaps") {
        std::string a = "ACGTACGT";
        std::string b = "TGCATGCA";
        int high_x = 10; // Increased mismatch penalty
        int expected_cost = -40; // 7I7D
        REQUIRE(wfa::naive(a, b, high_x, o, e) == expected_cost);
        REQUIRE(wfa::wavefront_dp(a, b, high_x, o, e) == expected_cost);
    }

    SECTION("Low Gap Penalty Preference for Mismatches") {
        std::string a = "AGCTAGCTAG";
        std::string b = "AGCTAG";
        int low_o = 1, low_e = 1; // Reduced gap penalties
        int expected_cost = -5; // 4D
        REQUIRE(wfa::naive(a, b, x, low_o, low_e) == expected_cost);
        REQUIRE(wfa::wavefront_dp(a, b, x, low_o, low_e) == expected_cost);
    }

    SECTION("Reversed Sequences") {
        std::string a = "GATTACA";
        std::string b = a;
        std::reverse(b.begin(), b.end());
        int expected_cost = -24; // 3X3X
        REQUIRE(wfa::naive(a, b, x, o, e) == expected_cost);
        REQUIRE(wfa::wavefront_dp(a, b, x, o, e) == expected_cost);
    }

    SECTION("Empty Sequence with Non-Zero Penalties") {
        std::string a = "";
        std::string b = "AGCT";
        int expected_cost = -14; // 4I
        REQUIRE(wfa::naive(a, b, x, o, e) == expected_cost);
        REQUIRE(wfa::wavefront_dp(a, b, x, o, e) == expected_cost);
    }

    SECTION("Sequences with Lowercase Letters") {
        std::string a = "gattaca";
        std::string b = "GATTACA";
        int expected_cost = -28; // 7X
        REQUIRE(wfa::naive(a, b, x, o, e) == expected_cost);
        REQUIRE(wfa::wavefront_dp(a, b, x, o, e) == expected_cost);
    }

    SECTION("Sequences with Special Characters") {
        std::string a = "GATTACA!";
        std::string b = "GATTACA?";
        int expected_cost = -4; // 1X
        REQUIRE(wfa::naive(a, b, x, o, e) == expected_cost);
        REQUIRE(wfa::wavefront_dp(a, b, x, o, e) == expected_cost);
    }

    SECTION("Very Long Sequences with Repeats") {
        std::string a(5000, 'A');
        std::string b(5000, 'A');
        a += std::string(5000, 'C');
        b += std::string(5000, 'C');
        int expected_cost = 0;
        REQUIRE(wfa::naive(a, b, x, o, e) == expected_cost);
        REQUIRE(wfa::wavefront_dp(a, b, x, o, e) == expected_cost);
    }
}

// Main function to run the Catch2 test session
int main(int argc, char* argv[]) {
    Catch::Session session; // Create a Catch2 session to run the tests

    // Run all test cases
    int result = session.run(argc, argv);

    // Return the result (ensure it's within valid range)
    return (result < 0xff ? result : 0xff);
}
