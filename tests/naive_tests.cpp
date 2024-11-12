#include "naive.hpp"          // Include the header for function declarations
#include <gtest/gtest.h>      // Google Test framework
#include <string>

// Fixture for Naive Alignment Tests
class NaiveAlignmentTest : public ::testing::Test {
protected:
    // Helper function to run the naive alignment and return the result
    int RunNaiveAlignment(const std::string& a, const std::string& b, int x, int o, int e) {
        return naive(a, b, x, o, e);
    }

    // Helper function to run the wavefront alignment and return the result
    int RunWavefrontAlignment(const std::string& a, const std::string& b, int x, int o, int e) {
        return wavefront(a, b, x, o, e);
    }
};

// Test Cases for Naive Alignment

// 1. Test identical sequences
TEST_F(NaiveAlignmentTest, IdenticalSequences) {
    std::string a = "GATTACA";
    std::string b = "GATTACA";
    int x = 4, o = 6, e = 2;

    EXPECT_EQ(RunNaiveAlignment(a, b, x, o, e), 0);
    EXPECT_EQ(RunWavefrontAlignment(a, b, x, o, e), 0);
}

// 2. Test completely different sequences
TEST_F(NaiveAlignmentTest, CompletelyDifferentSequences) {
    std::string a = "GATTACA";
    std::string b = "CTGACGT";
    int x = 4, o = 6, e = 2;

    EXPECT_GT(RunNaiveAlignment(a, b, x, o, e), 0);
    EXPECT_GT(RunWavefrontAlignment(a, b, x, o, e), 0);
}

// 3. Test single insertion penalty
TEST_F(NaiveAlignmentTest, SingleInsertion) {
    std::string a = "GATTACA";
    std::string b = "GATTTACA";
    int x = 4, o = 6, e = 2;

    int expected_cost = o + e;
    EXPECT_EQ(RunNaiveAlignment(a, b, x, o, e), expected_cost);
    EXPECT_EQ(RunWavefrontAlignment(a, b, x, o, e), expected_cost);
}

// 4. Test single deletion penalty
TEST_F(NaiveAlignmentTest, SingleDeletion) {
    std::string a = "GATTACA";
    std::string b = "GATA";
    int x = 4, o = 6, e = 2;

    int expected_cost = o + 3 * e;
    EXPECT_EQ(RunNaiveAlignment(a, b, x, o, e), expected_cost);
    EXPECT_EQ(RunWavefrontAlignment(a, b, x, o, e), expected_cost);
}

// 5. Test single mismatch penalty
TEST_F(NaiveAlignmentTest, SingleMismatch) {
    std::string a = "GATTACA";
    std::string b = "GACTACA";
    int x = 4, o = 6, e = 2;

    int expected_cost = x;
    EXPECT_EQ(RunNaiveAlignment(a, b, x, o, e), expected_cost);
    EXPECT_EQ(RunWavefrontAlignment(a, b, x, o, e), expected_cost);
}

// 6. Test longer sequences with minimal changes
TEST_F(NaiveAlignmentTest, LongSequenceWithSingleChange) {
    std::string a(1000, 'A');
    std::string b = a;
    b[500] = 'T';  // Introduce a single mismatch in the middle
    int x = 4, o = 6, e = 2;

    int expected_cost = x;
    EXPECT_EQ(RunNaiveAlignment(a, b, x, o, e), expected_cost);
    EXPECT_EQ(RunWavefrontAlignment(a, b, x, o, e), expected_cost);
}

// 7. Test empty sequences
TEST_F(NaiveAlignmentTest, EmptySequences) {
    std::string a = "";
    std::string b = "";
    int x = 4, o = 6, e = 2;

    EXPECT_EQ(RunNaiveAlignment(a, b, x, o, e), 0);
    EXPECT_EQ(RunWavefrontAlignment(a, b, x, o, e), 0);
}

// 8. Test one empty sequence
TEST_F(NaiveAlignmentTest, OneEmptySequence) {
    std::string a = "GATTACA";
    std::string b = "";
    int x = 4, o = 6, e = 2;

    int expected_cost = o + e * a.size();
    EXPECT_EQ(RunNaiveAlignment(a, b, x, o, e), expected_cost);
    EXPECT_EQ(RunWavefrontAlignment(a, b, x, o, e), expected_cost);
}

// 9. Test realistic DNA sequences with few differences
TEST_F(NaiveAlignmentTest, RealisticDNASequences) {
    std::string a = "AGCTTAGCTA";
    std::string b = "AGCTCGCTA";
    int x = 4, o = 6, e = 2;

    EXPECT_GT(RunNaiveAlignment(a, b, x, o, e), 0);
    EXPECT_GT(RunWavefrontAlignment(a, b, x, o, e), 0);
}

// 10. Test high penalties to check alignment behavior under extreme conditions
TEST_F(NaiveAlignmentTest, HighPenaltyConditions) {
    std::string a = "GATTACA";
    std::string b = "CTGACGT";
    int x = 10, o = 10, e = 10;

    EXPECT_GT(RunNaiveAlignment(a, b, x, o, e), 50);
    EXPECT_GT(RunWavefrontAlignment(a, b, x, o, e), 50);
}
