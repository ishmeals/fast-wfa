#pragma once

#include "wfa.hpp"

#include <string>
#include <string_view>
#include <vector>
#include <array>
#include <unordered_set>

#if _MSC_VER && !__INTEL_COMPILER && ! __clang__
#pragma warning(push)
#pragma warning(disable: 4996 4245 4324 4267 4244)
#include "Kokkos_SIMD.hpp"
#pragma warning(pop)
#else
#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wlanguage-extension-token"
#include "Kokkos_SIMD.hpp"
#pragma clang diagnostic pop
#endif

namespace wfa {
	//Native SIMD types will be optimised for the same architecture they are compiled on
	using simd_type = Kokkos::Experimental::native_simd<int32_t>;
	using mask_type = Kokkos::Experimental::native_simd_mask<int32_t>;
	using tag_type = Kokkos::Experimental::element_aligned_tag;

	//This function takes a wavefront entry, a column, a starting diagonal, and a scaling modifier to that diagonal (k + scaling), and will attempt to fill a vector register with the contents of the wavefront column entry in the range [k + scaling, k + scaling + simd_width). Out of bounds indicies are automatically handled, and replaced with -1 in the returned vector.
	simd_type simd_lookup(wavefront_entry_t& t, int32_t column, int32_t k, int32_t scaling);

	bool extend_simd(wavefront_t& wavefront, std::string_view a, std::string_view b);

	void next_simd(wavefront_t& wavefront, int32_t s, int32_t x, int32_t o, int32_t e);

	//The wavefront alignment algorithm implementation with manual vectorization. Aligns strings a and b with the provided substitution cost x, open cost o, and extend cost e. Uses the provided memory arena as it runs
	int32_t wavefront_simd(std::string_view a, std::string_view b, int32_t x, int32_t o, int32_t e, wavefront_arena_t& arena);


}