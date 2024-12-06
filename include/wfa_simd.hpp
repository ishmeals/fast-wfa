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
	using simd_type = Kokkos::Experimental::native_simd<int32_t>;
	using mask_type = Kokkos::Experimental::native_simd_mask<int32_t>;
	using tag_type = Kokkos::Experimental::element_aligned_tag;

	bool extend_simd(wavefront_t& wavefront, std::string_view a, std::string_view b);

	void next_simd(wavefront_t& wavefront, int32_t s, int32_t x, int32_t o, int32_t e);

	int32_t wavefront_simd(std::string_view a, std::string_view b, int32_t x, int32_t o, int32_t e);


}