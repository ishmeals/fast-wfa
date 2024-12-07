#include "include/data_gen.hpp"
#include "include/wfa_simd.hpp"
#include "include/wfa.hpp"
#include "include/naive.hpp"
#include "fmt/format.h"
#include "fmt/chrono.h"
#include <chrono>

//Naive is disabled for performance reasons

int main() {
    auto sequences = wfa::modify_sequences(100, 100000, 0.2);
    /*int32_t i = 0;
    wfa::wavefront_arena_t arena;
    for (const auto& pair : sequences) {
        const std::string& a = pair.first;
        const std::string& b = pair.second;
        auto score = wfa::wavefront_simd(a, b, 4, 6, 2, arena);
        auto score2 = wfa::naive(a, b, 4, 6, 2);
        auto score3 = wfa::wavefront(a, b, 4, 6, 2, arena);
        if (score != score2 || score != score3) {
            fmt::println("---Error---\nScore: {} {} {}\nA: {}\nB: {}", score, score2, score3, a, b);
        }
        else {
            fmt::println("Success: {}", i);
        }
        ++i;
    }*/

    auto start = std::chrono::system_clock::now();
    wfa::wavefront_arena_t arena1;
    for (const auto& pair : sequences) {
        const std::string& a = pair.first;
        const std::string& b = pair.second;
        wfa::wavefront(a, b, 4, 6, 2, arena1);
    }
    auto end = std::chrono::system_clock::now();
    fmt::println("Wavefront: {:%T}", end - start);

    start = std::chrono::system_clock::now();
    wfa::wavefront_arena_t arena2;
    for (const auto& pair : sequences) {
        const std::string& a = pair.first;
        const std::string& b = pair.second;
        wfa::wavefront_simd(a, b, 4, 6, 2, arena2);
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

}