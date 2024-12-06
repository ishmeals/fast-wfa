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

	struct wavefront_entry_t {
		int32_t low;
		int32_t high;
		int32_t number_per_col;
		std::vector<int32_t> data;
		bool valid;

		int32_t lookup(int32_t column, int32_t k);
		void set(int32_t column, int32_t k, int32_t val);
		int32_t& no_bound(int32_t column, int32_t row);
		int32_t* start_ptr(int32_t column);

		wavefront_entry_t(int32_t low, int32_t high);
		wavefront_entry_t();
	};

	struct wavefront_t {
		//std::vector<std::array<std::vector<int32_t>, 3>> data;
		std::vector<wavefront_entry_t> data;
		//std::vector<std::array<int32_t, 2>> low_hi;
		//ankerl::unordered_dense::set<int32_t> valid_scores;

		int32_t lookup(int32_t score, int32_t column, int32_t k);
		int32_t wave_size_low(int32_t score);
		int32_t wave_size_high(int32_t score);
		bool valid_score(int32_t score);
		
		void print();
	};

	//bounded score lookup
	//int32_t bn_s(wavefront_t& wavefront, int32_t score, int32_t column, int32_t row);

	bool extend(wavefront_t& wavefront, std::string_view a, std::string_view b);

	void next(wavefront_t& wavefront, int32_t s, int32_t x, int32_t o, int32_t e);

	int32_t wavefront(std::string_view a, std::string_view b, int32_t x, int32_t o, int32_t e);
}