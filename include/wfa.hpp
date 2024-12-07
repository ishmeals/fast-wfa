#pragma once

#include <string>
#include <string_view>
#include <vector>
#include <array>
#include <unordered_set>
#include "ankerl/unordered_dense.h"
#include <span>

namespace wfa {
	//using wavefront_t = std::vector<std::array<std::vector<int32_t>, 3>>;
	int32_t constexpr ins = 0;
	int32_t constexpr del = 1;
	int32_t constexpr match = 2;

	//A memory arena for use in the wave front types
	//Reuse between runs of the algorithm
	struct wavefront_arena_t {
		std::vector<int32_t> data;
		size_t current_index = 0;

		//Allocates <size> number of int32_t, returns a pointer to the beginning of the set
		int32_t* alloc(size_t size);

	};

	//An individual wavefront for the WFA algorithm
	struct wavefront_entry_t {
		//The range of diagonals covered
		int32_t low;
		int32_t high;
		int32_t number_per_col;
		
		//A view to the underlying arena allocation
		std::span<int32_t> data;

		int32_t lookup(int32_t column, int32_t k);
		void set(int32_t column, int32_t k, int32_t val);
		//Unbounded lookup - does not correct for row vs diagonal
		int32_t& no_bound(int32_t column, int32_t row);
		//returns the pointer to the start of the given column
		int32_t* start_ptr(int32_t column);

		wavefront_entry_t(int32_t low, int32_t high, wavefront_arena_t& arena);
		wavefront_entry_t(); //Default constructor for false states

		wavefront_entry_t(wavefront_entry_t&& rhs);
		wavefront_entry_t& operator=(wavefront_entry_t&& rhs);
		//Move only!
		wavefront_entry_t(const wavefront_entry_t& rhs) = delete;
		wavefront_entry_t& operator=(const wavefront_entry_t& rhs) = delete;
	};

	//The storage type for a full pass of the WFA algorithm
	struct wavefront_t {
		//A reference to the underlying memory arena
		wavefront_arena_t& arena;
		//A vector of wavefronts
		std::vector<wavefront_entry_t> views;
		//A vector mapping scores as indexs to indexs into views. -1 sentinal value for no wavefront with that score
		std::vector<int32_t> mapping;

		int32_t lookup(int32_t score, int32_t column, int32_t k);
		int32_t wave_size_low(int32_t score);
		int32_t wave_size_high(int32_t score);

		//checks if a score is valid and exists in the views
		bool valid_score(int32_t score);
		
		//inserts a new wavefront and handles reallocation and updating views
		wavefront_entry_t& insert(int32_t low, int32_t high);
		//inserts a dummy wavefront with no underlying allocation
		void insert();

		void print();

		wavefront_t(wavefront_arena_t& arena) : arena(arena) {}
	};

	bool extend(wavefront_t& wavefront, std::string_view a, std::string_view b);

	void next(wavefront_t& wavefront, int32_t s, int32_t x, int32_t o, int32_t e);


	//The wavefront alignment algorithm implementation. Aligns strings a and b with the provided substitution cost x, open cost o, and extend cost e. Uses the provided memory arena as it runs
	int32_t wavefront(std::string_view a, std::string_view b, int32_t x, int32_t o, int32_t e, wavefront_arena_t& arena);
}