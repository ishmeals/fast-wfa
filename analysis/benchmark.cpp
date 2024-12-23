#include "include/wfa.hpp"
#include "include/wfa_simd.hpp"
#include "include/naive.hpp"
#include "bindings/cpp/WFAligner.hpp"
#include "include/data_gen.hpp"
#include "fmt/format.h"
#include "fmt/chrono.h"
#include <chrono>
#include <fstream>
#include <cstdlib>
#include <string>

int main(int argc, char* argv[]) {
    // Check command-line arguments
    if (argc != 7) {
        fmt::println("Usage: {} <error_rate> <sequence_length> <num_sequences> <x> <o> <e>", argv[0]);
        return 1;
    }

    // Parse command-line arguments
    double error_rate = std::stod(argv[1]);
    int sequence_length = std::stoi(argv[2]);
    int num_sequences = std::stoi(argv[3]);
    int x = std::stoi(argv[4]);
    int o = std::stoi(argv[5]);
    int e = std::stoi(argv[6]);

    // Generate sequences
    auto sequences = wfa::modify_sequences(sequence_length, num_sequences, error_rate);

    // Benchmark Naive
    auto start = std::chrono::system_clock::now();
    wfa::wavefront_arena_t arena1;
    for (const auto& pair : sequences) {
        const std::string& a = pair.first;
        const std::string& b = pair.second;
        wfa::naive(a, b, x, o, e);
    }
    auto end = std::chrono::system_clock::now();
    fmt::println("Naive: {:%T}", end - start);


    // Benchmark Wavefront
    start = std::chrono::system_clock::now();
    for (const auto& pair : sequences) {
        const std::string& a = pair.first;
        const std::string& b = pair.second;
        wfa::wavefront(a, b, x, o, e, arena1);
    }
    end = std::chrono::system_clock::now();
    fmt::println("Wavefront: {:%T}", end - start);

    // Benchmark Wavefront SIMD
    start = std::chrono::system_clock::now();
    wfa::wavefront_arena_t arena2;
    for (const auto& pair : sequences) {
        const std::string& a = pair.first;
        const std::string& b = pair.second;
        wfa::wavefront_simd(a, b, x, o, e, arena2);
    }
    end = std::chrono::system_clock::now();
    fmt::println("Wavefront SIMD: {:%T}", end - start);

    // Benchmark WFA2-lib
    wfa::WFAlignerGapAffine aligner(x, o, e, wfa::WFAligner::Alignment, wfa::WFAligner::MemoryHigh);
    start = std::chrono::system_clock::now();
    for (const auto& pair : sequences) {
        const std::string& a = pair.first;
        const std::string& b = pair.second;

        aligner.alignEnd2End(a, b);
    }
    end = std::chrono::system_clock::now();
    fmt::println("WFA2-lib: {:%T}", end - start);
    return 0;
}
