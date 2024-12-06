#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <utility>
#include <cstdlib>
#include <sstream>
#include <chrono>
#include "include/data_gen.hpp"

// Function to execute the alignment executable and capture the output
std::tuple<double, double, double> run_alignment(const std::string& executable, const std::string& seq1, const std::string& seq2, int x, int o, int e) {
    std::ostringstream command;
    command << executable << " \"" << seq1 << "\" \"" << seq2 << "\" " << x << " " << o << " " << e;
    std::string cmd = command.str() + " > temp_output.txt";

    auto start = std::chrono::high_resolution_clock::now();
    int exit_code = std::system(cmd.c_str());
    auto end = std::chrono::high_resolution_clock::now();

    if (exit_code != 0) {
        std::cerr << "Command failed: " << cmd << "\n";
        return { -1.0, -1.0, -1.0 };
    }

    std::ifstream temp_output("temp_output.txt");
    if (!temp_output) {
        std::cerr << "Failed to open temp_output.txt\n";
        return { -1.0, -1.0, -1.0 };
    }

    std::string line;
    double avg_time = -1.0, avg_mem = -1.0, avg_utilization = -1.0;

    while (std::getline(temp_output, line)) {
        if (line.find("Execution Time:") != std::string::npos) {
            avg_time = std::stod(line.substr(line.find(":") + 1));
        }
        else if (line.find("Memory Consumption:") != std::string::npos) {
            avg_mem = std::stod(line.substr(line.find(":") + 1));
        }
        else if (line.find("Effective Core Utilization:") != std::string::npos) {
            avg_utilization = std::stod(line.substr(line.find(":") + 1));
        }
    }
    temp_output.close();
    return { avg_time, avg_mem, avg_utilization };
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

// A single experiment to vary sequence length
void experiment_vary_length(const std::string& executable, const std::string& output_csv) {
    int mismatch_penalty = 4;
    int gap_opening_cost = 6;
    int gap_extension_cost = 2;
    int num_samples = 100;
    double error_rate = 0.1;

    std::vector<int> sequence_lengths = { 100, 200, 500, 1000, 2000 };
    std::vector<std::vector<std::string>> results;

    for (int length : sequence_lengths) {
        auto sequences = wfa::modify_sequences(length, num_samples, error_rate);
        double total_time = 0.0, total_mem = 0.0, total_utilization = 0.0;

        for (const auto& [seq1, seq2] : sequences) {
            auto [time, mem, utilization] = run_alignment(executable, seq1, seq2, mismatch_penalty, gap_opening_cost, gap_extension_cost);
            total_time += time;
            total_mem += mem;
            total_utilization += utilization;
        }

        results.push_back({
            std::to_string(num_samples),
            std::to_string(length),
            std::to_string(error_rate),
            "0.5", // GC content assumed to be ~50% for random sequences
            std::to_string(mismatch_penalty),
            std::to_string(gap_opening_cost),
            std::to_string(gap_extension_cost),
            std::to_string(total_time / num_samples),
            std::to_string(total_mem / num_samples),
            std::to_string(total_utilization / num_samples)
            });
    }

    write_to_csv(output_csv, results);
}

// TODO: implement more experiments

// Main function
int main() {
    const std::string executable = "out/build/linux-debug/bin/wfa2_comparison";
    const std::string output_csv = "results/exp_results.csv";

    // Write CSV header
    std::ofstream csv_file(output_csv);
    csv_file << "Sample Count,Sequence Length,Error Rate,% GC Content,Mismatch Penalty,Gap Opening Cost,Gap Extension Cost,Avg Time,Avg Memory,Avg Core Utilization\n";
    csv_file.close();

    // Run experiments
    experiment_vary_length(executable, output_csv);

    // TODO: add additional experiment calls

    return 0;
}
