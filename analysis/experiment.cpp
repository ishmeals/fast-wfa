#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <cstdlib>
#include <sstream>
#include <chrono>
#include <iomanip>
#include <filesystem>
#include "include/data_gen.hpp"

// Helper function to extract time in seconds from the formatted output
double parse_time(const std::string& line) {
    if (line.empty()) return 0.0;
    size_t colon1 = line.find(':');
    size_t colon2 = line.find(':', colon1 + 1);
    size_t dot = line.find('.');

    if (colon1 == std::string::npos || colon2 == std::string::npos || dot == std::string::npos) {
        std::cerr << "Error parsing time: " << line << "\n";
        return 0.0;
    }

    int hours = std::stoi(line.substr(0, colon1));
    int minutes = std::stoi(line.substr(colon1 + 1, colon2 - colon1 - 1));
    double seconds = std::stod(line.substr(colon2 + 1));

    return hours * 3600 + minutes * 60 + seconds;
}

// Function to execute the alignment executable and capture the output
double run_alignment(const std::string& executable, double error_rate, int sequence_length,
    int num_sequences, int x, int o, int e, const std::string& algorithm) {
    std::ostringstream command;
    command << executable << " " << error_rate << " " << sequence_length << " " << num_sequences
        << " " << x << " " << o << " " << e;

    std::string cmd = command.str() + " > temp_output.txt";
    int exit_code = std::system(cmd.c_str());

    if (exit_code != 0) {
        std::cerr << "Command failed: " << cmd << "\n";
        return 0.0;
    }

    std::ifstream temp_output("temp_output.txt");
    if (!temp_output) {
        std::cerr << "Failed to open temp_output.txt\n";
        return 0.0;
    }

    std::string line;
    while (std::getline(temp_output, line)) {
        if (line.find(algorithm + ":") != std::string::npos) {
            return parse_time(line.substr(line.find(':') + 1));
        }
    }

    return 0.0; // Return 0 if time is not found
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

void experiment_vary_error_rate(const std::string& executable, const std::string& output_csv) {
    int mismatch_penalty = 4;
    int gap_opening_cost = 6;
    int gap_extension_cost = 2;
    int num_samples = 1000;
    int sequence_length = 100;
    std::vector<double> error_rates = { 0.01, 0.05, 0.1, 0.2, 0.3 };

    std::vector<std::string> algorithms = { "Wavefront", "Wavefront SIMD", "WFA2-lib" };
    for (double error_rate : error_rates) {
        for (const auto& algorithm : algorithms) {
            double avg_time = run_alignment(executable, error_rate, sequence_length, num_samples, mismatch_penalty, gap_opening_cost, gap_extension_cost, algorithm);
            write_to_csv(output_csv, { {algorithm, "Error v Time", std::to_string(num_samples), std::to_string(sequence_length), std::to_string(error_rate), std::to_string(mismatch_penalty), std::to_string(gap_opening_cost), std::to_string(gap_extension_cost), std::to_string(avg_time)} });
        }
    }
}

void experiment_vary_sequence_length(const std::string& executable, const std::string& output_csv) {
    double error_rate = 0.1;
    int mismatch_penalty = 4;
    int gap_opening_cost = 6;
    int gap_extension_cost = 2;
    int num_samples = 1000;
    std::vector<int> sequence_lengths = { 50, 100, 200, 500, 1000 };

    std::vector<std::string> algorithms = { "Wavefront", "Wavefront SIMD", "WFA2-lib" };
    for (int sequence_length : sequence_lengths) {
        for (const auto& algorithm : algorithms) {
            double avg_time = run_alignment(executable, error_rate, sequence_length, num_samples, mismatch_penalty, gap_opening_cost, gap_extension_cost, algorithm);
            write_to_csv(output_csv, { {algorithm, "Sequence Length v Time", std::to_string(num_samples), std::to_string(sequence_length), std::to_string(error_rate), std::to_string(mismatch_penalty), std::to_string(gap_opening_cost), std::to_string(gap_extension_cost), std::to_string(avg_time)} });
        }
    }
}

void experiment_vary_gap_opening(const std::string& executable, const std::string& output_csv) {
    double error_rate = 0.1;
    int sequence_length = 100;
    int mismatch_penalty = 4;
    int gap_extension_cost = 2;
    int num_samples = 1000;
    std::vector<int> gap_opening_costs = { 1, 3, 6, 10, 15 };

    std::vector<std::string> algorithms = { "Wavefront", "Wavefront SIMD", "WFA2-lib" };
    for (int gap_opening_cost : gap_opening_costs) {
        for (const auto& algorithm : algorithms) {
            double avg_time = run_alignment(executable, error_rate, sequence_length, num_samples, mismatch_penalty, gap_opening_cost, gap_extension_cost, algorithm);
            write_to_csv(output_csv, { {algorithm, "Gap Opening v Time", std::to_string(num_samples), std::to_string(sequence_length), std::to_string(error_rate), std::to_string(mismatch_penalty), std::to_string(gap_opening_cost), std::to_string(gap_extension_cost), std::to_string(avg_time)} });
        }
    }
}

void experiment_vary_gap_extension(const std::string& executable, const std::string& output_csv) {
    double error_rate = 0.1;
    int sequence_length = 100;
    int mismatch_penalty = 4;
    int gap_opening_cost = 6;
    int num_samples = 1000;
    std::vector<int> gap_extension_costs = { 1, 2, 5, 10, 20 };

    std::vector<std::string> algorithms = { "Wavefront", "Wavefront SIMD", "WFA2-lib" };
    for (int gap_extension_cost : gap_extension_costs) {
        for (const auto& algorithm : algorithms) {
            double avg_time = run_alignment(executable, error_rate, sequence_length, num_samples, mismatch_penalty, gap_opening_cost, gap_extension_cost, algorithm);
            write_to_csv(output_csv, { {algorithm, "Gap Extension v Time", std::to_string(num_samples), std::to_string(sequence_length), std::to_string(error_rate), std::to_string(mismatch_penalty), std::to_string(gap_opening_cost), std::to_string(gap_extension_cost), std::to_string(avg_time)} });
        }
    }
}

void experiment_vary_mismatch_penalty(const std::string& executable, const std::string& output_csv) {
    double error_rate = 0.1;
    int sequence_length = 100;
    int gap_opening_cost = 6;
    int gap_extension_cost = 2;
    int num_samples = 1000;
    std::vector<int> mismatch_penalties = { 1, 2, 4, 8, 10 };

    std::vector<std::string> algorithms = { "Wavefront", "Wavefront SIMD", "WFA2-lib" };
    for (int mismatch_penalty : mismatch_penalties) {
        for (const auto& algorithm : algorithms) {
            double avg_time = run_alignment(executable, error_rate, sequence_length, num_samples, mismatch_penalty, gap_opening_cost, gap_extension_cost, algorithm);
            write_to_csv(output_csv, { {algorithm, "Mismatch Penalty v Time", std::to_string(num_samples), std::to_string(sequence_length), std::to_string(error_rate), std::to_string(mismatch_penalty), std::to_string(gap_opening_cost), std::to_string(gap_extension_cost), std::to_string(avg_time)} });
        }
    }
}

// Experiment: Joint Impact of Error Rate and Sequence Length
void experiment_joint_error_length(const std::string& executable, const std::string& output_csv) {
    int mismatch_penalty = 4;
    int gap_opening_cost = 6;
    int gap_extension_cost = 2;
    int num_samples = 1000;
    std::vector<int> sequence_lengths = { 50, 100, 200, 500 };
    std::vector<double> error_rates = { 0.01, 0.05, 0.1, 0.2 };

    std::vector<std::string> algorithms = { "Wavefront", "Wavefront SIMD", "WFA2-lib" };
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
void experiment_interaction_gap_costs(const std::string& executable, const std::string& output_csv) {
    double error_rate = 0.1;
    int sequence_length = 100;
    int mismatch_penalty = 4;
    int num_samples = 1000;
    std::vector<int> gap_opening_costs = { 1, 3, 6, 10 };
    std::vector<int> gap_extension_costs = { 1, 2, 5, 10 };

    std::vector<std::string> algorithms = { "Wavefront", "Wavefront SIMD", "WFA2-lib" };
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
void experiment_sensitivity_analysis(const std::string& executable, const std::string& output_csv) {
    double error_rate = 0.1;
    int sequence_length = 100;
    int num_samples = 1000;
    std::vector<int> mismatch_penalties = { 2, 4, 8 };
    std::vector<int> gap_opening_costs = { 3, 6, 10 };
    std::vector<int> gap_extension_costs = { 1, 5, 10 };

    std::vector<std::string> algorithms = { "Wavefront", "Wavefront SIMD", "WFA2-lib" };
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
void experiment_error_rate_complexity(const std::string& executable, const std::string& output_csv) {
    int sequence_length = 100;
    int mismatch_penalty = 4;
    int gap_extension_cost = 2;
    int num_samples = 1000;
    std::vector<double> error_rates = { 0.01, 0.05, 0.1, 0.2, 0.3 };
    std::vector<int> gap_opening_costs = { 6, 10 };

    std::vector<std::string> algorithms = { "Wavefront", "Wavefront SIMD", "WFA2-lib" };
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
void experiment_length_gap_penalties(const std::string& executable, const std::string& output_csv) {
    double error_rate = 0.05;
    int mismatch_penalty = 4;
    int gap_extension_cost = 2;
    int num_samples = 1000;
    std::vector<int> sequence_lengths = { 50, 100, 200, 500 };
    std::vector<int> gap_opening_costs = { 3, 6, 10 };

    std::vector<std::string> algorithms = { "Wavefront", "Wavefront SIMD", "WFA2-lib" };
    for (int sequence_length : sequence_lengths) {
        for (int gap_opening_cost : gap_opening_costs) {
            for (const auto& algorithm : algorithms) {
                double avg_time = run_alignment(executable, error_rate, sequence_length, num_samples, mismatch_penalty, gap_opening_cost, gap_extension_cost, algorithm);
                write_to_csv(output_csv, { {algorithm, "Length & Gap Penalties", std::to_string(num_samples), std::to_string(sequence_length), std::to_string(error_rate), std::to_string(mismatch_penalty), std::to_string(gap_opening_cost), std::to_string(gap_extension_cost), std::to_string(avg_time)} });
            }
        }
    }
}

// Helper function to find the repository base
std::string get_repository_base() {
    std::filesystem::path current_path = std::filesystem::current_path();
    while (!current_path.empty() && !std::filesystem::exists(current_path / ".git")) {
        current_path = current_path.parent_path();
    }
    if (current_path.empty()) {
        throw std::runtime_error("Repository base could not be determined. Make sure you are running within a git repository.");
    }
    return current_path.string();
}

// Main function
int main() {
    // Dynamically determine the repository base and paths
    std::string repo_base = get_repository_base();
    std::string executable = repo_base + "/out/build/linux-debug/bin/wfa2_comparison";
    std::string output_csv = repo_base + "/results/exp_results.csv";

    // Ensure the results directory exists
    std::filesystem::create_directories(repo_base + "/results");
    const std::string executable = "~/fast-wfa/out/build/linux-debug/bin/wfa2_comparison";
    const std::string output_csv = "~/fast-wfa/results/exp_results.csv";

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
        experiment_vary_sequence_length(executable, output_csv);

        std::cout << "Running Experiment: Gap Opening Cost vs Time...\n";
        experiment_vary_gap_opening(executable, output_csv);

        std::cout << "Running Experiment: Gap Extension Cost vs Time...\n";
        experiment_vary_gap_extension(executable, output_csv);

        std::cout << "Running Experiment: Mismatch Penalty vs Time...\n";
        experiment_vary_mismatch_penalty(executable, output_csv);

        std::cout << "Running Experiment: Joint Impact of Error Rate and Sequence Length...\n";
        experiment_joint_error_length(executable, output_csv);

        std::cout << "Running Experiment: Interaction of Gap Costs...\n";
        experiment_interaction_gap_costs(executable, output_csv);

        std::cout << "Running Experiment: Sensitivity Analysis...\n";
        experiment_sensitivity_analysis(executable, output_csv);

        std::cout << "Running Experiment: Varying Error Rate with Increased Complexity...\n";
        experiment_error_rate_complexity(executable, output_csv);

        std::cout << "Running Experiment: Combination of Sequence Length and Gap Penalties...\n";
        experiment_length_gap_penalties(executable, output_csv);

        std::cout << "All experiments completed. Results written to: " << output_csv << "\n";
    }
    catch (const std::exception& e) {
        std::cerr << "An error occurred during execution: " << e.what() << "\n";
        return 1;
    }

    return 0;
}

