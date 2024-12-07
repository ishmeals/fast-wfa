# WFA

This is an implementation of the Wavefront Alignment algorithm in both a autovectorized approach and a manually vectorized approach. This repository also contains the functionality to compare the implementations, alongside comparing with the venerable WFA2-lib.

# Building

Given the fragile nature of building C++ project, we provide a dockerfile that is capable of building and runnning the program as it is. While this is the recommended approach, manual build steps are listed below.

## CMake

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

`wfa_tool` is setup to run a large number of sequences, and automatically print the results. 

`wfa2_comparison` can be run with 

