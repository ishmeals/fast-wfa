#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <chrono>
#include <iomanip>
#include "include/wfa.hpp"
#include "include/wfa_simd.hpp"
#include "include/naive.hpp"
#include "bindings/cpp/WFAligner.hpp"
#include "include/data_gen.hpp"

// Function to execute and benchmark algorithms
double run_alignment(double error_rate, int sequence_length, int num_sequences,
    int mismatch_penalty, int gap_opening_cost, int gap_extension_cost,
    const std::string& algorithm) {
    // Generate sequences
    auto sequences = wfa::modify_sequences(sequence_length, num_sequences, error_rate);

    auto start = std::chrono::high_resolution_clock::now();
    if (algorithm == "Naive") {
        for (const auto& pair : sequences) {
            const std::string& a = pair.first;
            const std::string& b = pair.second;
            wfa::naive(a, b, mismatch_penalty, gap_opening_cost, gap_extension_cost);
        }
    }
    else if (algorithm == "Wavefront") {
        wfa::wavefront_arena_t arena;
        for (const auto& pair : sequences) {
            const std::string& a = pair.first;
            const std::string& b = pair.second;
            wfa::wavefront(a, b, mismatch_penalty, gap_opening_cost, gap_extension_cost, arena);
        }
    }
    else if (algorithm == "Wavefront SIMD") {
        wfa::wavefront_arena_t arena;
        for (const auto& pair : sequences) {
            const std::string& a = pair.first;
            const std::string& b = pair.second;
            wfa::wavefront_simd(a, b, mismatch_penalty, gap_opening_cost, gap_extension_cost, arena);
        }
    }
    else if (algorithm == "WFA2-lib") {
        wfa::WFAlignerGapAffine aligner(mismatch_penalty, gap_opening_cost, gap_extension_cost,
            wfa::WFAligner::Alignment, wfa::WFAligner::MemoryHigh);
        for (const auto& pair : sequences) {
            const std::string& a = pair.first;
            const std::string& b = pair.second;
            aligner.alignEnd2End(a, b);
        }
    }
    auto end = std::chrono::high_resolution_clock::now();

    // Calculate and return elapsed time in milliseconds
    return std::chrono::duration<double, std::milli>(end - start).count();
}

// Function to write data to the CSV file
void write_to_csv(const std::string& filename, const std::vector<std::vector<std::string>>& data) {
    std::ofstream csv_file(filename, std::ios::app);
    if (!csv_file) {
        std::cerr << "Failed to open file: " << filename << "\n";
        return;
    }

    for (const auto& row : data) {
        for (size_t i = 0; i < row.size(); ++i) {
            csv_file << row[i];
            if (i < row.size() - 1) csv_file << ",";
        }
        csv_file << "\n";
    }

    csv_file.close();
}


void experiment_vary_error_rate(const std::string& output_csv) {
    int mismatch_penalty = 4;
    int gap_opening_cost = 6;
    int gap_extension_cost = 2;
    int num_samples = 1000;
    int sequence_length = 100;
    std::vector<double> error_rates = { 0.01, 0.05, 0.1, 0.2, 0.3 };

    std::vector<std::string> algorithms = { "Naive", "Wavefront", "Wavefront SIMD", "WFA2-lib" };
    for (double error_rate : error_rates) {
        for (const auto& algorithm : algorithms) {
            double avg_time = run_alignment(executable, error_rate, sequence_length, num_samples, mismatch_penalty, gap_opening_cost, gap_extension_cost, algorithm);
            write_to_csv(output_csv, { {algorithm, "Error v Time", std::to_string(num_samples), std::to_string(sequence_length), std::to_string(error_rate), std::to_string(mismatch_penalty), std::to_string(gap_opening_cost), std::to_string(gap_extension_cost), std::to_string(avg_time)} });
        }
    }
}

void experiment_vary_sequence_length(const std::string& output_csv) {
    double error_rate = 0.1;
    int mismatch_penalty = 4;
    int gap_opening_cost = 6;
    int gap_extension_cost = 2;
    int num_samples = 1000;
    std::vector<int> sequence_lengths = { 50, 100, 200, 500, 1000 };

    std::vector<std::string> algorithms = { "Naive", "Wavefront", "Wavefront SIMD", "WFA2-lib" };
    for (int sequence_length : sequence_lengths) {
        for (const auto& algorithm : algorithms) {
            double avg_time = run_alignment(executable, error_rate, sequence_length, num_samples, mismatch_penalty, gap_opening_cost, gap_extension_cost, algorithm);
            write_to_csv(output_csv, { {algorithm, "Sequence Length v Time", std::to_string(num_samples), std::to_string(sequence_length), std::to_string(error_rate), std::to_string(mismatch_penalty), std::to_string(gap_opening_cost), std::to_string(gap_extension_cost), std::to_string(avg_time)} });
        }
    }
}

void experiment_vary_gap_opening(const std::string& output_csv) {
    double error_rate = 0.1;
    int sequence_length = 100;
    int mismatch_penalty = 4;
    int gap_extension_cost = 2;
    int num_samples = 1000;
    std::vector<int> gap_opening_costs = { 1, 3, 6, 10, 15 };

    std::vector<std::string> algorithms = { "Naive", "Wavefront", "Wavefront SIMD", "WFA2-lib" };
    for (int gap_opening_cost : gap_opening_costs) {
        for (const auto& algorithm : algorithms) {
            double avg_time = run_alignment(executable, error_rate, sequence_length, num_samples, mismatch_penalty, gap_opening_cost, gap_extension_cost, algorithm);
            write_to_csv(output_csv, { {algorithm, "Gap Opening v Time", std::to_string(num_samples), std::to_string(sequence_length), std::to_string(error_rate), std::to_string(mismatch_penalty), std::to_string(gap_opening_cost), std::to_string(gap_extension_cost), std::to_string(avg_time)} });
        }
    }
}

void experiment_vary_gap_extension(const std::string& output_csv) {
    double error_rate = 0.1;
    int sequence_length = 100;
    int mismatch_penalty = 4;
    int gap_opening_cost = 6;
    int num_samples = 1000;
    std::vector<int> gap_extension_costs = { 1, 2, 5, 10, 20 };

    std::vector<std::string> algorithms = { "Naive", "Wavefront", "Wavefront SIMD", "WFA2-lib" };
    for (int gap_extension_cost : gap_extension_costs) {
        for (const auto& algorithm : algorithms) {
            double avg_time = run_alignment(executable, error_rate, sequence_length, num_samples, mismatch_penalty, gap_opening_cost, gap_extension_cost, algorithm);
            write_to_csv(output_csv, { {algorithm, "Gap Extension v Time", std::to_string(num_samples), std::to_string(sequence_length), std::to_string(error_rate), std::to_string(mismatch_penalty), std::to_string(gap_opening_cost), std::to_string(gap_extension_cost), std::to_string(avg_time)} });
        }
    }
}

void experiment_vary_mismatch_penalty(const std::string& output_csv) {
    double error_rate = 0.1;
    int sequence_length = 100;
    int gap_opening_cost = 6;
    int gap_extension_cost = 2;
    int num_samples = 1000;
    std::vector<int> mismatch_penalties = { 1, 2, 4, 8, 10 };

    std::vector<std::string> algorithms = { "Naive", "Wavefront", "Wavefront SIMD", "WFA2-lib" };
    for (int mismatch_penalty : mismatch_penalties) {
        for (const auto& algorithm : algorithms) {
            double avg_time = run_alignment(executable, error_rate, sequence_length, num_samples, mismatch_penalty, gap_opening_cost, gap_extension_cost, algorithm);
            write_to_csv(output_csv, { {algorithm, "Mismatch Penalty v Time", std::to_string(num_samples), std::to_string(sequence_length), std::to_string(error_rate), std::to_string(mismatch_penalty), std::to_string(gap_opening_cost), std::to_string(gap_extension_cost), std::to_string(avg_time)} });
        }
    }
}

// Experiment: Joint Impact of Error Rate and Sequence Length
void experiment_joint_error_length(const std::string& output_csv) {
    int mismatch_penalty = 4;
    int gap_opening_cost = 6;
    int gap_extension_cost = 2;
    int num_samples = 1000;
    std::vector<int> sequence_lengths = { 50, 100, 200, 500 };
    std::vector<double> error_rates = { 0.01, 0.05, 0.1, 0.2 };

    std::vector<std::string> algorithms = { "Naive", "Wavefront", "Wavefront SIMD", "WFA2-lib" };
    for (int sequence_length : sequence_lengths) {
        for (double error_rate : error_rates) {
            for (const auto& algorithm : algorithms) {
                double avg_time = run_alignment(executable, error_rate, sequence_length, num_samples, mismatch_penalty, gap_opening_cost, gap_extension_cost, algorithm);
                write_to_csv(output_csv, { {algorithm, "Joint Error & Length", std::to_string(num_samples), std::to_string(sequence_length), std::to_string(error_rate), std::to_string(mismatch_penalty), std::to_string(gap_opening_cost), std::to_string(gap_extension_cost), std::to_string(avg_time)} });
            }
        }
    }
}

// Experiment: Interaction of Gap Costs
void experiment_interaction_gap_costs(const std::string& output_csv) {
    double error_rate = 0.1;
    int sequence_length = 100;
    int mismatch_penalty = 4;
    int num_samples = 1000000;
    std::vector<int> gap_opening_costs = { 1, 3, 6, 10 };
    std::vector<int> gap_extension_costs = { 1, 2, 5, 10 };

    std::vector<std::string> algorithms = { "Naive", "Wavefront", "Wavefront SIMD", "WFA2-lib" };
    for (int gap_opening_cost : gap_opening_costs) {
        for (int gap_extension_cost : gap_extension_costs) {
            for (const auto& algorithm : algorithms) {
                double avg_time = run_alignment(executable, error_rate, sequence_length, num_samples, mismatch_penalty, gap_opening_cost, gap_extension_cost, algorithm);
                write_to_csv(output_csv, { {algorithm, "Gap Costs Interaction", std::to_string(num_samples), std::to_string(sequence_length), std::to_string(error_rate), std::to_string(mismatch_penalty), std::to_string(gap_opening_cost), std::to_string(gap_extension_cost), std::to_string(avg_time)} });
            }
        }
    }
}

// Experiment: Algorithm Sensitivity Analysis
void experiment_sensitivity_analysis(const std::string& output_csv) {
    double error_rate = 0.1;
    int sequence_length = 100;
    int num_samples = 1000000;
    std::vector<int> mismatch_penalties = { 2, 4, 8 };
    std::vector<int> gap_opening_costs = { 3, 6, 10 };
    std::vector<int> gap_extension_costs = { 1, 5, 10 };

    std::vector<std::string> algorithms = { "Naive", "Wavefront", "Wavefront SIMD", "WFA2-lib" };
    for (int mismatch_penalty : mismatch_penalties) {
        for (int gap_opening_cost : gap_opening_costs) {
            for (int gap_extension_cost : gap_extension_costs) {
                for (const auto& algorithm : algorithms) {
                    double avg_time = run_alignment(executable, error_rate, sequence_length, num_samples, mismatch_penalty, gap_opening_cost, gap_extension_cost, algorithm);
                    write_to_csv(output_csv, { {algorithm, "Sensitivity Analysis", std::to_string(num_samples), std::to_string(sequence_length), std::to_string(error_rate), std::to_string(mismatch_penalty), std::to_string(gap_opening_cost), std::to_string(gap_extension_cost), std::to_string(avg_time)} });
                }
            }
        }
    }
}

// Experiment: Varying Error Rate with Increased Complexity
void experiment_error_rate_complexity(const std::string& output_csv) {
    int sequence_length = 100;
    int mismatch_penalty = 4;
    int gap_extension_cost = 2;
    int num_samples = 1000;
    std::vector<double> error_rates = { 0.01, 0.05, 0.1, 0.2, 0.3 };
    std::vector<int> gap_opening_costs = { 6, 10 };

    std::vector<std::string> algorithms = { "Naive", "Wavefront", "Wavefront SIMD", "WFA2-lib" };
    for (double error_rate : error_rates) {
        for (int gap_opening_cost : gap_opening_costs) {
            for (const auto& algorithm : algorithms) {
                double avg_time = run_alignment(executable, error_rate, sequence_length, num_samples, mismatch_penalty, gap_opening_cost, gap_extension_cost, algorithm);
                write_to_csv(output_csv, { {algorithm, "Error Rate & Complexity", std::to_string(num_samples), std::to_string(sequence_length), std::to_string(error_rate), std::to_string(mismatch_penalty), std::to_string(gap_opening_cost), std::to_string(gap_extension_cost), std::to_string(avg_time)} });
            }
        }
    }
}

// Experiment: Combination of Sequence Length and Gap Penalties
void experiment_length_gap_penalties(const std::string& output_csv) {
    double error_rate = 0.05;
    int mismatch_penalty = 4;
    int gap_extension_cost = 2;
    int num_samples = 1000;
    std::vector<int> sequence_lengths = { 50, 100, 200, 500 };
    std::vector<int> gap_opening_costs = { 3, 6, 10 };

    std::vector<std::string> algorithms = { "Naive", "Wavefront", "Wavefront SIMD", "WFA2-lib" };
    for (int sequence_length : sequence_lengths) {
        for (int gap_opening_cost : gap_opening_costs) {
            for (const auto& algorithm : algorithms) {
                double avg_time = run_alignment(executable, error_rate, sequence_length, num_samples, mismatch_penalty, gap_opening_cost, gap_extension_cost, algorithm);
                write_to_csv(output_csv, { {algorithm, "Length & Gap Penalties", std::to_string(num_samples), std::to_string(sequence_length), std::to_string(error_rate), std::to_string(mismatch_penalty), std::to_string(gap_opening_cost), std::to_string(gap_extension_cost), std::to_string(avg_time)} });
            }
        }
    }
}


// Main function
int main() {
    std::string output_csv = "exp_results.csv";

    // Write CSV header
    std::ofstream csv_file(output_csv);
    if (!csv_file.is_open()) {
        std::cerr << "Failed to open file: " << output_csv << "\n";
        return 1;
    }

    csv_file << "Algorithm,Experiment,Sample Count,Sequence Length,Error Rate,Mismatch Penalty,Gap Opening Cost,Gap Extension Cost,Avg Time\n";
    csv_file.close();

    try {
        // Run all experiments
        std::cout << "Running Experiment: Error Rate vs Time...\n";
        experiment_vary_error_rate(executable, output_csv);

        std::cout << "Running Experiment: Sequence Length vs Time...\n";
        // experiment_vary_sequence_length(executable, output_csv);

        std::cout << "Running Experiment: Gap Opening Cost vs Time...\n";
        // experiment_vary_gap_opening(executable, output_csv);

        std::cout << "Running Experiment: Gap Extension Cost vs Time...\n";
        // experiment_vary_gap_extension(executable, output_csv);

        std::cout << "Running Experiment: Mismatch Penalty vs Time...\n";
        // experiment_vary_mismatch_penalty(executable, output_csv);

        std::cout << "Running Experiment: Joint Impact of Error Rate and Sequence Length...\n";
        // experiment_joint_error_length(executable, output_csv);

        std::cout << "Running Experiment: Interaction of Gap Costs...\n";
        // experiment_interaction_gap_costs(executable, output_csv);

        std::cout << "Running Experiment: Sensitivity Analysis...\n";
        // experiment_sensitivity_analysis(executable, output_csv);

        std::cout << "Running Experiment: Varying Error Rate with Increased Complexity...\n";
        experiment_error_rate_complexity(executable, output_csv);

        std::cout << "Running Experiment: Combination of Sequence Length and Gap Penalties...\n";
        // experiment_length_gap_penalties(executable, output_csv);

        std::cout << "All experiments completed. Results written to: " << output_csv << "\n";
    }
    catch (const std::exception& e) {
        std::cerr << "An error occurred during execution: " << e.what() << "\n";
        return 1;
    }

    return 0;
}