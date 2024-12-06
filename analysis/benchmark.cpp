#include "include/wfa.hpp"
#include "include/wfa_simd.hpp"
#include "bindings/cpp/WFAligner.hpp"
#include "include/data_gen.hpp"
#include "fmt/format.h"
#include "fmt/chrono.h"
#include <chrono>

int main() {
	
	auto sequences = wfa::modify_sequences(100, 100000, 0.01);

    auto start = std::chrono::system_clock::now();
    wfa::wavefront_arena_t arena1;
    for (const auto& pair : sequences) {
        const std::string& a = pair.first;
        const std::string& b = pair.second;
        wfa::wavefront(a, b, 4, 6, 2, arena1);
        //fmt::println("--");
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

    wfa::WFAlignerGapAffine aligner(4, 6, 2, wfa::WFAligner::Alignment, wfa::WFAligner::MemoryHigh);
    start = std::chrono::system_clock::now();
    for (const auto& pair : sequences) {
        const std::string& a = pair.first;
        const std::string& b = pair.second;
        
        aligner.alignEnd2End(a, b);
    }
    end = std::chrono::system_clock::now();
    fmt::println("WFA2-lib: {:%T}", end - start);
}