#include "include/wfa.hpp"
#include "include/wfa_simd.hpp"
#include "bindings/cpp/WFAligner.hpp"
#include "include/data_gen.hpp"
#include "fmt/format.h"
#include "fmt/chrono.h"
#include <chrono>

int main() {
	
	auto sequences = wfa::modify_sequences(100, 10000, 0.05);

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

    start = std::chrono::system_clock::now();
    for (const auto& pair : sequences) {
        const std::string& a = pair.first;
        const std::string& b = pair.second;
        wfa::WFAlignerGapAffine aligner(4, 6, 2, wfa::WFAligner::Alignment, wfa::WFAligner::MemoryHigh);
        aligner.alignEnd2End(a, b);
    }
    end = std::chrono::system_clock::now();
    fmt::println("WFA2-lib: {:%T}", end - start);
}