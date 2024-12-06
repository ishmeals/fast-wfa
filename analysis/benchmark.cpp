#include <iostream>
#include <string>
#include <string_view>
#include <chrono>
#include <vector>
#include <array>
#include <algorithm>

#include "include/naive.hpp"
#include "include/wfa.hpp"
#include "include/wfa_simd.hpp"
#include "bindings/cpp/WFAligner.hpp"

wavefront_aligner_t* init_aligner(int x, int o, int e) {
    wavefront_aligner_attr_t attributes = wavefront_aligner_attr_default;
    attributes.distance_metric = gap_affine;
    attributes.affine_penalties.mismatch = x;
    attributes.affine_penalties.gap_opening = o;
    attributes.affine_penalties.gap_extension = e;
    attributes.alignment_scope = compute_score;

    return wavefront_aligner_new(&attributes);
}

int wfalib2_align(std::string_view a, std::string_view b, wavefront_aligner_t* wf_aligner) {
    std::string pattern(a);
    std::string text(b);

    wavefront_aligner_t* wf_aligner = wavefront_aligner_new(&attributes);
    wavefront_align(wf_aligner, pattern.c_str(), (int)pattern.size(), text.c_str(), (int)text.size());
    return wf_aligner->cigar.score;
}

void run_vtune(const std::string& algorithm, const std::string& output_file, const std::string& command) {
    std::string vtune_command = "vtune -collect hotspots -result-dir vtune_results/" + algorithm + " -- " + command;
    std::cout << "Running VTune for algorithm: " << algorithm << "\n";
    std::system(vtune_command.c_str());
    std::cout << "VTune results saved for " << algorithm << " in vtune_results/" + algorithm + "\n";
}

// Benchmarking function
template <typename Func>
void benchmark_algorithm(const std::string& algorithm_name, Func align_func, const std::string& seq1, const std::string& seq2, int x, int o, int e) {
    auto start = std::chrono::high_resolution_clock::now();

    // Run the alignment algorithm
    align_func(seq1, seq2, x, o, e);

    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed_time = end - start;

    std::cout << "Algorithm: " << algorithm_name << "\n";
    std::cout << "Execution Time: " << elapsed_time.count() << " seconds\n";

    std::string command = "./" + algorithm_name + "_binary";
    run_vtune(algorithm_name, algorithm_name + "_results", command);
}

int main(int argc, char* argv[]) {
    if (argc < 6) {
        std::cerr << "Usage: " << argv[0] << " <seq1> <seq2> <x> <o> <e> [alg]\n";
        return 1;
    }

    std::string seq1 = argv[1];
    std::string seq2 = argv[2];
    int x = std::atoi(argv[3]);
    int o = std::atoi(argv[4]);
    int e = std::atoi(argv[5]);

    wavefront_aligner_t* aligner = init_aligner(x, o, e);

    if (argc == 7) {
        std::string selected_algorithm = argv[6];
        if (selected_algorithm == "naive") benchmark_algorithm("Naive", wfa::naive, seq1, seq2, x, o, e);
        else if (selected_algorithm == "dp_wfa") benchmark_algorithm("DP_WFA", wfa::wavefront_dp, seq1, seq2, x, o, e);
        else if (selected_algorithm == "trad_wfa") benchmark_algorithm("Traditional_WFA", wfa::wavefront, seq1, seq2, x, o, e);
        // else if (selected_algorithm == "simd_wfa") benchmark_algorithm("SIMD_WFA", simd_wfa_align, seq1, seq2, x, o, e);
        else if (selected_algorithm == "wfa2_lib") benchmark_algorithm("WFA2_Lib", wfalib2_align, seq1, seq2, x, o, e);
        else {
            std::cerr << "Unknown algorithm specified.\n";
            return 1;
        }
    }
    else {
        // Run all algorithms
        benchmark_algorithm("Naive", wfa::naive, seq1, seq2, x, o, e);
        benchmark_algorithm("DP_WFA", wfa::wavefront_dp, seq1, seq2, x, o, e);
        benchmark_algorithm("Traditional_WFA", wfa::wavefront, seq1, seq2, x, o, e);
        // benchmark_algorithm("SIMD_WFA", simd_wfa_align, seq1, seq2, x, o, e);
        benchmark_algorithm("WFA2_Lib", wfalib2_align, seq1, seq2, x, o, e);
    }

    return 0;
}