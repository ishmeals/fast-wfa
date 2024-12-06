#include "include/data_gen.hpp"
#include "include/wfa_simd.hpp"
#include "include/wfa.hpp"
#include "include/naive.hpp"
#include "fmt/format.h"
#include "fmt/chrono.h"
#include <chrono>

int main() {
    auto sequences = wfa::modify_sequences(100, 100000, 0.02);
    /*int32_t i = 0;
    for (const auto& pair : sequences) {
        const std::string& a = pair.first;
        const std::string& b = pair.second;
        auto score = wfa::wavefront_simd(a, b, 4, 6, 2);
        auto score2 = wfa::naive(a, b, 4, 6, 2);
        auto score3 = wfa::wavefront(a, b, 4, 6, 2);
        if (score != score2 || score != score3) {
            fmt::println("---Error---\nScore: {} {} {}\nA: {}\nB: {}", score, score2, score3, a, b);
        }
        else {
            fmt::println("Sucess: {}", i);
        }
        ++i;
    }*/

    auto start = std::chrono::system_clock::now();
    for (const auto& pair : sequences) {
        const std::string& a = pair.first;
        const std::string& b = pair.second;
        wfa::wavefront(a, b, 4, 6, 2);
    }
    auto end = std::chrono::system_clock::now();
    fmt::println("Wavefront: {:%T}", end - start);

    start = std::chrono::system_clock::now();
    for (const auto& pair : sequences) {
        const std::string& a = pair.first;
        const std::string& b = pair.second;
        wfa::wavefront_simd(a, b, 4, 6, 2);
    }
    end = std::chrono::system_clock::now();
    fmt::println("Wavefront SIMD: {:%T}", end - start);

    /*start = std::chrono::system_clock::now();
    for (const auto& pair : sequences) {
        const std::string& a = pair.first;
        const std::string& b = pair.second;
        wfa::naive(a, b, 4, 6, 2);
    }
    end = std::chrono::system_clock::now();
    fmt::println("Naive: {:%T}", end - start);*/

   /* start = std::chrono::system_clock::now();
    for (const auto& pair : sequences) {
        const std::string& a = pair.first;
        const std::string& b = pair.second;
        wfa::wavefront(a, b, 4, 6, 2);
    }
    end = std::chrono::system_clock::now();
    fmt::println("Wavefront DP: {:%T}", end - start);*/


    /*std::string a = "TCATGCAGGCCCGGAGGAAAATTATCCCGGGACGATCTCTAGTTATTGCAGTGGCCAGCGTTGGGGACGGTCGTTTTCGCTCCCAGAACGCCAACCGTGT";
    std::string b = "TCATGCAGGCCCGAAGGAAAATTAGTCCGGGATGATTTCTTGCTATTTCAGTGGCCAGCATTGGGGCATATCGCTTCAGTTCCCATAAGGTCGACCGTGA";

    auto score = wfa::wavefront_simd(a, b, 4, 6, 2);
    auto score2 = wfa::naive(a, b, 4, 6, 2);
    auto score3 = wfa::wavefront(a, b, 4, 6, 2);

    fmt::println("{} {} {}", score, score2, score3);*/
   /*  std::string a = "AGGATGCTCG";
    std::string b = "ACCATACTCG";
    
    fmt::println("{} {} {}", score, score2, score3);*/
}