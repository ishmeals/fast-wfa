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

// Function to align using WFA2-lib
void wfalib2_align(const std::string& seq1, const std::string& seq2, int x, int o, int e) {
    wfa::WFAlignerGapAffine aligner(x, o, e, wfa::WFAligner::Alignment, wfa::WFAligner::MemoryHigh);
    aligner.alignEnd2End(seq1, seq2);
    std::cout << "Alignment score (WFA2): " << aligner.getAlignmentScore() << "\n";
}

// Function to run VTune profiler
void run_vtune(const std::string& algorithm, const std::string& executable,
    const std::string& seq1, const std::string& seq2, int x, int o, int e) {
    std::string vtune_dir = "results/" + algorithm;
    try {
        std::filesystem::create_directories(vtune_dir);
    }
    catch (const std::filesystem::filesystem_error& e) {
        std::cerr << "Failed to create VTune results directory: " << e.what() << "\n";
        return;
    }

    std::string vtune_command = "vtune -collect hotspots -result-dir " + vtune_dir +
        " -- " + executable + " " + seq1 + " " + seq2 + " " +
        std::to_string(x) + " " + std::to_string(o) + " " +
        std::to_string(e) + " " + algorithm;
    std::cout << "Running VTune for algorithm: " << algorithm << "\n";
    int result = std::system(vtune_command.c_str());
    if (result != 0) {
        std::cerr << "VTune command failed with exit code: " << result << "\n";
    }
    else {
        std::cout << "VTune results saved in: " << vtune_dir << "\n";
    }
}

// Function to benchmark a single algorithm
template <typename Func>
void benchmark_algorithm(const std::string& algorithm_name, Func align_func,
    const std::string& seq1, const std::string& seq2, int x, int o, int e) {
    auto start = std::chrono::high_resolution_clock::now();
    align_func(seq1, seq2, x, o, e);
    auto end = std::chrono::high_resolution_clock::now();

    std::chrono::duration<double> elapsed_time = end - start;
    std::cout << "Algorithm: " << algorithm_name << "\n";
    std::cout << "Execution Time: " << elapsed_time.count() << " seconds\n";
}

int main(int argc, char* argv[]) {
    if (argc < 6) {
        std::cerr << "Usage: " << argv[0] << " <seq1> <seq2> <x> <o> <e> [algorithm]\n";
        return 1;
    }

    // Parse inputs
    std::string seq1 = argv[1];
    std::string seq2 = argv[2];
    int x = std::stoi(argv[3]);
    int o = std::stoi(argv[4]);
    int e = std::stoi(argv[5]);
    std::string executable = argv[0];

    // First, run all algorithms if no algorithm is specified
    if (argc == 6) {
        benchmark_algorithm("Naive", wfa::naive, seq1, seq2, x, o, e);
        benchmark_algorithm("DP_WFA", wfa::wavefront_dp, seq1, seq2, x, o, e);
        benchmark_algorithm("Traditional_WFA", wfa::wavefront, seq1, seq2, x, o, e);
        benchmark_algorithm("WFA2_Lib", wfalib2_align, seq1, seq2, x, o, e);
    }
    // Then, run VTune on a specific algorithm if specified
    else if (argc == 7) {
        std::string selected_algorithm = argv[6];
        if (selected_algorithm == "naive") {
            run_vtune("Naive", executable, seq1, seq2, x, o, e);
        }
        else if (selected_algorithm == "dp_wfa") {
            run_vtune("DP_WFA", executable, seq1, seq2, x, o, e);
        }
        else if (selected_algorithm == "trad_wfa") {
            run_vtune("Traditional_WFA", executable, seq1, seq2, x, o, e);
        }
        else if (selected_algorithm == "wfa2_lib") {
            run_vtune("WFA2_Lib", executable, seq1, seq2, x, o, e);
        }
        else {
            std::cerr << "Unknown algorithm: " << selected_algorithm << "\n";
            return 1;
        }
    }

    return 0;
}
