#pragma once

#include "kokkos_simd.hpp"

#include <string>
#include <string_view>
#include <vector>
#include <array>
#include <unordered_map>

namespace wfa {
	//using wavefront_t = std::vector<std::array<std::vector<int32_t>, 3>>;
	int32_t constexpr match = 2;
	int32_t constexpr ins = 0;
	int32_t constexpr del = 1;

	struct wavefront_t {
		std::vector<std::array<std::vector<int32_t>, 3>> data;
		std::unordered_map<int32_t, int32_t> score_to_index;
		std::vector<std::array<int32_t, 2>> low_hi;

		int32_t lookup(int32_t score, int32_t column, int32_t row);
		int32_t wave_size(int32_t score, bool low);

	};

	//bounded score lookup
	//int32_t bn_s(wavefront_t& wavefront, int32_t score, int32_t column, int32_t row);

	bool extend(wavefront_t& wavefront, std::string_view a, std::string_view b, int32_t score);

	void next(wavefront_t& wavefront, int32_t s, int32_t x, int32_t o, int32_t e);

	int32_t wavefront_simd(std::string_view a, std::string_view b, int32_t x, int32_t o, int32_t e);


}