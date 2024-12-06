#pragma once

#include <string>
#include <string_view>
#include <vector>
#include <array>
#include <unordered_set>
#include "ankerl/unordered_dense.h"

namespace wfa {
	//using wavefront_t = std::vector<std::array<std::vector<int32_t>, 3>>;
	int32_t constexpr ins = 0;
	int32_t constexpr del = 1;
	int32_t constexpr match = 2;

	struct wavefront_t {
		std::vector<std::array<std::vector<int32_t>, 3>> data;
		std::vector<std::array<int32_t, 2>> low_hi;
		ankerl::unordered_dense::set<int32_t> valid_scores;

		int32_t lookup(int32_t score, int32_t column, int32_t k);
		int32_t wave_size(int32_t score, bool low);

		void print();
	};

	//bounded score lookup
	//int32_t bn_s(wavefront_t& wavefront, int32_t score, int32_t column, int32_t row);

	bool extend(wavefront_t& wavefront, std::string_view a, std::string_view b, int32_t score);

	void next(wavefront_t& wavefront, int32_t s, int32_t x, int32_t o, int32_t e);

	int32_t wavefront(std::string_view a, std::string_view b, int32_t x, int32_t o, int32_t e);
}