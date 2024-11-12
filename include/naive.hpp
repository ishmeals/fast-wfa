#ifndef NAIVE_HPP
#define NAIVE_HPP

#include <string>
#include <string_view>

/**
 * @brief Calculates the alignment cost between two sequences using a naive dynamic programming approach.
 * 
 * The `naive` function computes the optimal alignment cost by filling a 2D matrix
 * that represents possible alignments between two input sequences. This method
 * is straightforward but can be inefficient for large sequences due to its O(n*m) complexity.
 * 
 * @param a The first sequence as a `std::string_view`.
 * @param b The second sequence as a `std::string_view`.
 * @param x The penalty cost for mismatches between sequences.
 * @param o The cost for opening a gap in the alignment.
 * @param e The cost for extending a gap in the alignment.
 * @return int The alignment cost between the two sequences.
 */
int naive(std::string_view a, std::string_view b, int x, int o, int e);

/**
 * @brief Computes the alignment cost between two sequences using a wavefront alignment approach.
 * 
 * The `wavefront` function optimizes the alignment calculation by focusing on the diagonals in the
 * alignment matrix, reducing unnecessary computations. This approach is generally more efficient
 * for larger sequences and may provide better performance than the naive method.
 * 
 * @param a The first sequence as a `std::string_view`.
 * @param b The second sequence as a `std::string_view`.
 * @param x The penalty cost for mismatches between sequences.
 * @param o The cost for opening a gap in the alignment.
 * @param e The cost for extending a gap in the alignment.
 * @return int The alignment cost between the two sequences.
 */
int wavefront(std::string_view a, std::string_view b, int x, int o, int e);

#endif // NAIVE_HPP
