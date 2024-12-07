# WFA

This is an implementation of the Wavefront Alignment algorithm in both a autovectorized approach and a manually vectorized approach. This repository also contains the functionality to compare the implementations, alongside comparing with the venerable WFA2-lib.

## Features
- Multiple alignment algorithms: Naive, Wavefront, SIMD Wavefront, and WFA2-lib integration.
- Experimentation framework for analyzing algorithmic performance under various conditions.
- Benchmarking tools for performance evaluation.
- Visualization scripts to generate performance graphs.
- Dockerized build and runtime environment.

# Building

Given the fragile nature of building C++ project, we provide a dockerfile that is capable of building and runnning the program as it is. While this is the recommended approach, manual build steps are listed below.

## Dockerfile
The dockerfile is setup to display benchmarking results depending based on the passed in parameters.

From the root of the repository, build the dockerfile with

```
docker build -t wfa .
```
Once the image is built, you can run 

```
docker run wfa <error_rate> <sequence_length> <num_sequences> <x> <o> <e>
```
For example, the command
```
docker run wfa 0.02 100 10000 4 6 2
```
Yields the following output
```
Naive: 00:00:00.130498
Wavefront: 00:00:00.002120
Wavefront SIMD: 00:00:00.002345
WFA2-lib: 00:00:00.003775
```

## Manual Approach

### Dependencies

This project requires the following dependencies:
- C++20 compiler (Clang 19+ recommended)
- CMake 3.28 or higher
- Ninja build system
- Python 3.8+ with pandas, seaborn, matplotlib
- Libraries (automatically pulled through FetchContent):
  - [Kokkos](https://github.com/kokkos/kokkos)
  - [fmt](https://github.com/fmtlib/fmt)
  - [WFA2-lib](https://github.com/QuantumFelidae/WFA2-lib)
  - [Catch2](https://github.com/catchorg/Catch2)



### Linux

A current version of clang-19 is the expected compiler on linux. If you do not have clang-19, it can be installed through the following series of commands

```
wget https://apt.llvm.org/llvm.sh
chmod +x llvm.sh
sudo ./llvm.sh 19 all
```

This project also uses CMake and Ninja as the build system. If your package manager does not have CMake 3.28 available, then current versions of CMake and Ninja can be retrieved through pip with the command:

```
python3 -m pip install --upgrade pip
python3 -m pip install cmake
python3 -m pip install ninja 
```

This project uses CMakePresets to automatically handle build parameters. There are two presets available on linux: `linux-debug` and `linux-release`. From the root directory of the repository, run

```
cmake --preset=linux-release
```

to build the release version of the repository. Build files will be placed in `./out/build/linux-release` or `./out/build/linux-debug`, depending on the preset.

To build all executables, enter the above build folder, and run

```
ninja
```

The executable files will be placed in `./out/build/<PRESET_NAME>/bin`

### Windows

If you have Visual Studio (NOT CODE) installed, then you have effectively everything needed to build the project. Note that you won't have access to WFA2-lib files, as there is no compatibility on windows.

The available presets are `x64-debug` and `x64-release`. If you installed clang-tools when you installed Visual Studio, then you will also have access to `x64-debug-clang` and `x64-release-clang`, which are recommended over the MSVC implementation for performance reasons.

# Code Structure

## File Structure

The `analysis` folder contains the C++ and python files used to generate the plots and data for the final paper.

The `include` and `src` files contain the main body of our WFA2 libary implementation.

The `tests` folder contains the set of tests used to check our implementation.

## Code structure

All of our library code is contained within the `wfa` namespace. The `naive.h/cpp` files contain the implementation of the SWG approach, and a naive dynamic programming approach to the WFA algorithm used for testing purposes. `wfa.h/cpp` contains the implementation of the wavefront data structures, and the basic `wfa::wavefront` implementation. `wfa_simd.h/cpp` contains the implementation of `wfa::simd__wavefront`. `data_gen.h/cpp` contains the code to generate the synthetic sequences.

# Executables

## Quick Tools

`wfa_tool` is setup to run a large number of random sequences, and automatically print the results. 

`wfa2_comparison` can be run with custom parameters to do benchmarks, see docker section above for details:

## Experiment
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
