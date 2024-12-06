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

// Experiment to vary error rate and record times
void experiment_vary_error_rate(const std::string& executable, const std::string& output_csv) {
    int mismatch_penalty = 4;
    int gap_opening_cost = 6;
    int gap_extension_cost = 2;
    int num_samples = 1000;
    int sequence_length = 100;

    std::vector<double> error_rates = { 0.01, 0.05, 0.1, 0.2, 0.3 };
    std::vector<std::string> algorithms = { "Wavefront", "Wavefront SIMD", "WFA2-lib" };
    std::vector<std::vector<std::string>> results;

    for (double error_rate : error_rates) {
        for (const auto& algorithm : algorithms) {
            double avg_time = run_alignment(executable, error_rate, sequence_length, num_samples,
                mismatch_penalty, gap_opening_cost, gap_extension_cost, algorithm);

            results.push_back({
                algorithm,
                "Error v Time",
                std::to_string(num_samples),
                std::to_string(sequence_length),
                std::to_string(error_rate),
                std::to_string(mismatch_penalty),
                std::to_string(gap_opening_cost),
                std::to_string(gap_extension_cost),
                std::to_string(avg_time)
                });
        }
    }

    write_to_csv(output_csv, results);
}

// Main function
int main() {
    std::filesystem::create_directories(std::filesystem::path("~/fast-wfa/results"));
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

    // Run experiment
    experiment_vary_error_rate(executable, output_csv);

    std::cout << "Experiment completed. Results written to: " << output_csv << "\n";

    return 0;
}
