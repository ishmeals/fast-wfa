#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <cstdlib>
#include <sstream>
#include <unistd.h>
#include <sys/types.h>
#include <pwd.h>

// Function to execute the alignment executable and capture the output time
double run_alignment(const std::string& executable, double error_rate, int sequence_length, int num_sequences, int x, int o, int e) {
    std::ostringstream command;
    command << executable << " " << error_rate << " " << sequence_length << " " << num_sequences << " " << x << " " << o << " " << e;
    std::string cmd = command.str() + " > temp_output.txt";

    int exit_code = std::system(cmd.c_str());
    if (exit_code != 0) {
        std::cerr << "Command failed: " << cmd << "\n";
        return -1.0;
    }

    std::ifstream temp_output("temp_output.txt");
    if (!temp_output) {
        std::cerr << "Failed to open temp_output.txt\n";
        return -1.0;
    }

    std::string line;
    double avg_time = -1.0;

    while (std::getline(temp_output, line)) {
        if (line.find("Wavefront:") != std::string::npos) {
            avg_time = std::stod(line.substr(line.find(":") + 1));
        }
    }

    temp_output.close();
    return avg_time;
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

// Experiment to vary error rate
void experiment_vary_error_rate(const std::string& executable, const std::string& output_csv) {
    int mismatch_penalty = 4;
    int gap_opening_cost = 6;
    int gap_extension_cost = 2;
    int num_sequences = 1000;
    int sequence_length = 100;

    std::vector<double> error_rates = { 0.01, 0.05, 0.1, 0.2, 0.3 };
    std::vector<std::vector<std::string>> results;

    for (double error_rate : error_rates) {
        double avg_time = run_alignment(executable, error_rate, sequence_length, num_sequences, mismatch_penalty, gap_opening_cost, gap_extension_cost);
        if (avg_time < 0) {
            std::cerr << "Error during alignment at error rate " << error_rate << "\n";
            continue;
        }

        results.push_back({
            std::to_string(num_sequences),
            std::to_string(sequence_length),
            std::to_string(error_rate),
            std::to_string(mismatch_penalty),
            std::to_string(gap_opening_cost),
            std::to_string(gap_extension_cost),
            std::to_string(avg_time)
            });
    }

    write_to_csv(output_csv, results);
}

// Main function
int main() {
    const std::string executable = "~/fast-wfa/out/build/linux-debug/bin/wfa2_comparison";
    const std::string output_csv = "exp_results.csv";

    // Write CSV header
    std::ofstream csv_file(output_csv);
    csv_file << "Sample Count,Sequence Length,Error Rate,Mismatch Penalty,Gap Opening Cost,Gap Extension Cost,Avg Time\n";
    csv_file.close();

    // Run experiment
    experiment_vary_error_rate(executable, output_csv);

    return 0;
}
