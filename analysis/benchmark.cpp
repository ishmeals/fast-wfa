#include <iostream>
#include <string>
#include <string_view>
#include <chrono>
#include <cstdlib>
#include <filesystem>

#include "include/naive.hpp"
#include "include/wfa.hpp"
#include "include/wfa_simd.hpp"
#include "bindings/cpp/WFAligner.hpp"

// Align using WFA2-lib
void wfalib2_align(const std::string& seq1, const std::string& seq2, int x, int o, int e) {
    wfa::WFAlignerGapAffine aligner(x, o, e, wfa::WFAligner::Alignment, wfa::WFAligner::MemoryHigh);
    aligner.alignEnd2End(seq1, seq2);
}

// Run VTune profiler
void run_vtune(const std::string& algorithm, const std::string& command) {
    std::string vtune_dir = "results/" + algorithm;
    std::filesystem::create_directories(vtune_dir); // Ensure directory exists
    std::string vtune_command = "vtune -collect hotspots -result-dir " + vtune_dir + " -- " + command;

    std::cout << "Running VTune for algorithm: " << algorithm << "\n";
    int result = std::system(vtune_command.c_str());
    if (result != 0) {
        std::cerr << "VTune command failed with exit code: " << result << "\n";
    }
    else {
        std::cout << "VTune results saved in: " << vtune_dir << "\n";
    }
}

// Benchmarking function
template <typename Func>
void benchmark_algorithm(const std::string& algorithm_name, Func align_func,
    const std::string& seq1, const std::string& seq2, int x, int o, int e) {
    // Measure time
    auto start = std::chrono::high_resolution_clock::now();
    align_func(seq1, seq2, x, o, e);
    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed_time = end - start;

    // Output runtime
    std::cout << "Algorithm: " << algorithm_name << "\n";
    std::cout << "Execution Time: " << elapsed_time.count() << " seconds\n";

    // VTune profiling
    std::string command = "./wfa2_comparison " + seq1 + " " + seq2 + " " + std::to_string(x)
        + " " + std::to_string(o) + " " + std::to_string(e) + " " + algorithm_name;
    run_vtune(algorithm_name, command);
}

int main(int argc, char* argv[]) {
    if (argc < 7) {
        std::cerr << "Usage: " << argv[0] << " <seq1> <seq2> <x> <o> <e> [algorithm]\n";
        return 1;
    }

    // Parse inputs
    std::string seq1 = argv[1];
    std::string seq2 = argv[2];
    int x = std::stoi(argv[3]);
    int o = std::stoi(argv[4]);
    int e = std::stoi(argv[5]);

    // Vtune running
    if (argc == 7) {
        std::string selected_algorithm = argv[6];
        if (selected_algorithm == "naive") {
            wfa::naive(seq1, seq2, x, o, e);
        }
        else if (selected_algorithm == "dp_wfa") {
            wfa::wavefront_dp(seq1, seq2, x, o, e);
        }
        else if (selected_algorithm == "trad_wfa") {
            wfa::wavefront(seq1, seq2, x, o, e);
        }
        else if (selected_algorithm == "wfa2_lib") {
            wfalib2_align(seq1, seq2, x, o, e);
        }

    }
    else {
        // Benchmarks
        benchmark_algorithm("Naive", wfa::naive, seq1, seq2, x, o, e);
        benchmark_algorithm("DP_WFA", wfa::wavefront_dp, seq1, seq2, x, o, e);
        benchmark_algorithm("Traditional_WFA", wfa::wavefront, seq1, seq2, x, o, e);
        benchmark_algorithm("WFA2_Lib", wfalib2_align, seq1, seq2, x, o, e);
    }

    return 0;
}
