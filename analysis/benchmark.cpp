#include <iostream>
#include <string>
#include <string_view>
#include <chrono>
#include <cstdlib>
#include <filesystem>
#include <fstream>

#include "include/naive.hpp"
#include "include/wfa.hpp"
#include "include/wfa_simd.hpp"
#include "bindings/cpp/WFAligner.hpp"

// Function to align using WFA2-lib
void wfalib2_align(const std::string& seq1, const std::string& seq2, int x, int o, int e) {
    wfa::WFAlignerGapAffine aligner(x, o, e, wfa::WFAligner::Alignment, wfa::WFAligner::MemoryHigh);
    aligner.alignEnd2End(seq1, seq2);
}

// Function to benchmark and time an algorithm
template <typename Func>
void benchmark_algorithm(const std::string& algorithm_name, Func align_func,
    const std::string& seq1, const std::string& seq2, int x, int o, int e) {
    // Measure runtime
    auto start = std::chrono::high_resolution_clock::now();
    align_func(seq1, seq2, x, o, e);
    auto end = std::chrono::high_resolution_clock::now();

    std::chrono::duration<double> elapsed_time = end - start;
    std::cout << "Algorithm: " << algorithm_name << "\n";
    std::cout << "Execution Time: " << elapsed_time.count() << " seconds\n";
}

// Function to invoke VTune profiling
void run_vtune(const std::string& executable, const std::string& seq1, const std::string& seq2,
    int x, int o, int e, const std::string& algorithm) {
    std::string vtune_dir = "results/" + algorithm;
    try {
        std::filesystem::create_directories(vtune_dir);
    }
    catch (const std::filesystem::filesystem_error& e) {
        std::cerr << "Failed to create VTune results directory: " << e.what() << "\n";
        return;
    }

    std::string vtune_command = "vtune -collect hotspots --summary -- " + executable + " " + seq1 + " " + seq2 + " " +
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

int main(int argc, char* argv[]) {

    std::ofstream output_file("output.txt");
    if (!output_file.is_open()) {
        std::cerr << "Failed to open output file.\n";
        return 1;
    }

    // Redirect std::cout and std::cerr to the file
    std::streambuf* cout_buf = std::cout.rdbuf(output_file.rdbuf());
    std::streambuf* cerr_buf = std::cerr.rdbuf(output_file.rdbuf());

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

    // If algorithm is specified, run it for VTune profiling only
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
        else {
            std::cerr << "Unknown algorithm: " << selected_algorithm << "\n";
            return 1;
        }
        return 0;
    }

    // If no algorithm is specified, benchmark all and invoke VTune profiling
    benchmark_algorithm("Naive", wfa::naive, seq1, seq2, x, o, e);
    run_vtune(executable, seq1, seq2, x, o, e, "naive");

    benchmark_algorithm("DP_WFA", wfa::wavefront_dp, seq1, seq2, x, o, e);
    run_vtune(executable, seq1, seq2, x, o, e, "dp_wfa");

    benchmark_algorithm("Traditional_WFA", wfa::wavefront, seq1, seq2, x, o, e);
    run_vtune(executable, seq1, seq2, x, o, e, "trad_wfa");

    benchmark_algorithm("WFA2_Lib", wfalib2_align, seq1, seq2, x, o, e);
    run_vtune(executable, seq1, seq2, x, o, e, "wfa2_lib");


    std::cout.rdbuf(cout_buf);
    std::cerr.rdbuf(cerr_buf);

    return 0;
}
