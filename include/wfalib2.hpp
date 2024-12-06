#pragma once
#include <string_view>

/**
 * @brief Perform gap-affine alignment using WFA2-lib.
 *
 * This function aligns two sequences using WFA2-lib's gap-affine mode.
 * @param a First sequence as a std::string_view
 * @param b Second sequence as a std::string_view
 * @param x Mismatch penalty (must be > 0)
 * @param o Gap-opening penalty (must be >= 0)
 * @param e Gap-extension penalty (must be > 0)
 * @return int The alignment score computed by WFA2-lib
 */
	int wfalib2_align(std::string_view a, std::string_view b, int x, int o, int e);
