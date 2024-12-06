#pragma once

#include "wfa.hpp"
#include "wfa_simd.hpp"

#include <string>
#include <string_view>
#include <vector>
#include <array>
#include <thread>
#include <semaphore>

#if _MSC_VER && !__INTEL_COMPILER && ! __clang__
#pragma warning(push)
#pragma warning(disable: 4996 4245 4324 4267 4244)
#include "Kokkos_Core.hpp"
#pragma warning(pop)
#else
#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wlanguage-extension-token"
#include "Kokkos_Core.hpp"
#pragma clang diagnostic pop
#endif

namespace wfa {
	struct resources_t {
		std::vector<wavefront_arena_t> arenas;
		std::vector<std::binary_semaphore> owners;
	};

	std::vector<int32_t> parallel(std::shared_ptr<std::vector<std::pair<std::string, std::string>>> sequences, int32_t x, int32_t o, int32_t e, bool use_simd, int32_t resources);
}
