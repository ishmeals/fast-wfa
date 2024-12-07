# Fast-WFA: High-Performance Sequence Alignment Library

Fast-WFA is a high-performance library for pairwise sequence alignment using the Wavefront Alignment Algorithm (WFA). This project includes several implementations of alignment algorithms, benchmarks, and tools for experimentation. It leverages SIMD optimization, hardware accelerations, and modern libraries like Kokkos for portability and performance.

Noah Edmiston, Logan Miller, Ismeal Lee
---

## Table of Contents
- [Features](#features)
- [Dependencies](#dependencies)
- [Installation](#installation)
- [Usage](#usage)
  - [Running Experiments](#running-experiments)
  - [Benchmarking](#benchmarking)
  - [Testing](#testing)
- [Directory Structure](#directory-structure)
- [Authors](#authors)
- [License](#license)

---

## Features

- Multiple alignment algorithms: Naive, Wavefront, SIMD Wavefront, and WFA2-lib integration.
- Experimentation framework for analyzing algorithmic performance under various conditions.
- Benchmarking tools for performance evaluation.
- Visualization scripts to generate performance graphs.
- Dockerized build and runtime environment.

---

## Dependencies

This project requires the following dependencies:
- C++20 compiler (Clang 19+ recommended)
- CMake 3.28 or higher
- Ninja build system
- Python 3.8+ with pandas, seaborn, matplotlib
- Libraries:
  - [Kokkos](https://github.com/kokkos/kokkos)
  - [fmt](https://github.com/fmtlib/fmt)
  - [WFA2-lib](https://github.com/QuantumFelidae/WFA2-lib)
  - [Catch2](https://github.com/catchorg/Catch2)

---

## Installation

### Building Locally
1. Clone the repository:
   ```bash
   git clone https://github.com/yourusername/fast-wfa.git
   cd fast-wfa
   ```
2. Build using CMake:
   ```bash
   cmake --preset=linux-release
   cmake --build ./out/build/linux-release
   ```


Executables:
### Benchmarking Tool
- Benchmarks the time performance for the following sequence alignment algorithms:
  - **Naive**
  - **Wavefront**
  - **Wavefront SIMD**
  - **WFA2-lib**
- Supports multiple experiments, such as varying:
  - Error rate
  - Sequence length
  - Gap opening and extension penalties
  - Mismatch penalties
  - Algorithmic complexity
- Outputs run parameters and results to a CSV file

#### How to Use

1. **Build the Project**
   - Create a build directory and compile the project:
     ```bash
     cmake --preset=linux-release
     cd out/build/linux-release
     ninja
     ```
   - Ensure you have all required dependencies installed (e.g., `fmt`, `WFA2-lib`).

2. **Run the Benchmarking Tool**
   - Execute the benchmarking tool to run predefined experiments:
     ```bash
     ./bin/experiment
     ```
   - By default, this will:
     - Vary parameters such as error rates, sequence lengths, and penalties.
     - Test the performance of multiple algorithms.
     - Generate a CSV file named `exp_results.csv` in the directory where the tool is executed.

3. **Interpret Results**
   - The tool generates a CSV file in the following format:
     ```csv
     Algorithm,Experiment,Sample Count,Sequence Length,Error Rate,Mismatch Penalty,Gap Opening Cost,Gap Extension Cost,Avg Time
     ```
   - This file can be used to analyze the performance of the algorithms under various conditions.

4. **Visualize Results**
   - Use the `graph.py` script to generate visualizations from the CSV file. For example:
     ```bash
     python graph.py /path/to/exp_results.csv
     ```
   - This script creates line plots and heatmaps based on the experimental data.

### Testing
1. Clone the repository:
   ```bash
   cmake --build ./out/build/linux-release --target naive_tests
   ./out/build/linux-release/bin/wfa_tests
   ```


## Directory Structure

.
├── CMakeLists.txt          # Project-wide CMake configuration
├── Dockerfile              # Dockerfile for building and running the project
├── include/                # Header files
│   ├── data_gen.hpp        # Sequence generation functions
│   ├── naive.hpp           # Naive DP algorithm
│   ├── wfa.hpp             # Wavefront Alignment algorithm
│   ├── wfa_simd.hpp        # SIMD-optimized Wavefront Alignment
├── src/                    # Source code
│   ├── naive.cpp           # Implementation of naive DP
│   ├── wfa.cpp             # Implementation of WFA
│   ├── wfa_simd.cpp        # SIMD implementation of WFA
│   ├── data_gen.cpp        # Sequence generation utilities
├── analysis/               # Experimentation and visualization
│   ├── experiment.cpp      # Experiment framework
│   ├── benchmark.cpp       # Benchmarking tool
│   ├── graph.py            # Visualization script
├── tests/                  # Unit tests
│   ├── naive_tests.cpp     # Catch2 tests for algorithms
│   ├── CMakeLists.txt      # Test build configuration
├── .git/                   # Git repository metadata
└── out/                    # Build output directory