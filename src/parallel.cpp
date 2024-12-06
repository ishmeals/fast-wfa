#include "include/parallel.hpp"

std::vector<int32_t> wfa::parallel(std::shared_ptr<std::vector<std::pair<std::string, std::string>>> sequences, int32_t x, int32_t o, int32_t e, bool use_simd, int32_t resources) {
	std::shared_ptr r = std::make_shared<resources_t>();
	r->arenas.resize(resources);
	for (int32_t i = 0; i < resources; ++i) {
		r->owners.push_back(std::binary_semaphore(1));
	}
	//r->owners.resize(resources, std::binary_semaphore(1));
	std::shared_ptr v = std::make_shared<std::vector<int32_t>>(sequences->size());

	auto f = KOKKOS_LAMBDA(size_t i) {
		auto pair = sequences->at(i);
		const std::string& a = pair.first;
		const std::string& b = pair.second;
		size_t r_i = 0;
		bool owned = false;
		while (!owned) {
			for (; r_i < r->owners.size(); ++i) {
				if (r->owners.at(r_i).try_acquire()) {
					owned = true;
					break;
				}
			}
		}

		v->at(i) = wfa::wavefront(a, b, x, o, e, r->arenas.at(r_i));
		r->owners.at(r_i).release();
	};

	Kokkos::parallel_for(sequences->size(), f);
}