#include "include/kokkos_simd.hpp"

#include "fmt/format.h"
#include "fmt/ranges.h"

int32_t wfa::wavefront_t::lookup(int32_t score, int32_t column, int32_t row) {
	if (score < 0) {
		return 0;
	}
	if (row < 0) {
		return 0;
	}
	auto iter = score_to_index.find(score);
	if (iter == score_to_index.end()) {
		return 0;
	}
	if (row > static_cast<int32_t>(data[iter->second][column].size()) - 1) {
		return 0;
	}
	return (data[iter->second][column][row]);
}

int32_t wfa::wavefront_t::wave_size(int32_t score, bool low) {
	auto iter = score_to_index.find(score);
	if (iter == score_to_index.end()) {
		return 0;
	}
	auto& pair = low_hi[iter->second];
	if (low) {
		return pair[0];
	}
	return pair[1];
}

void wfa::wavefront_t::print() {
	for (const auto& pair : score_to_index) {
		fmt::println("Score: {}", pair.first);
		fmt::println("I: {}\nD: {}\nM: {}", data[pair.second][ins], data[pair.second][del], data[pair.second][match]);
	}
}
//int32_t wfa::bn_s(wavefront_t& wavefront, int32_t score, int32_t column, int32_t row) {
//	if (score < 0) {
//		return 0;
//	}
//	return wavefront[score][column][row];
//}

bool wfa::extend(wavefront_t& wavefront, std::string_view a, std::string_view b, int32_t score) {
	std::vector<int32_t>& matchfront_back = wavefront.data.back()[2];
	int32_t k_low = wavefront.wave_size(score, true);
	int32_t k_high = wavefront.wave_size(score, false);
	for (int32_t k = k_low; k < k_high + 1; ++k) {
		int32_t starting_index = matchfront_back[k - k_low];
		int32_t v = starting_index - k;
		int32_t h = starting_index;
		if (v >= static_cast<int32_t>(a.size()) or h >= static_cast<int32_t>(b.size())) {
			continue;
		}
		bool mismatch = false;
		while (not mismatch) {
			const char v_c = a.at(v);
			const char h_c = b.at(h);
			if (v_c == h_c) {
				if (v == static_cast<int32_t>(a.size()) - 1 and h == static_cast<int32_t>(b.size()) - 1) {
					return true;
				}
				++starting_index;
				++v;
				++h;
			}
			else {
				mismatch = true;
			}
		}
		matchfront_back[k - k_low] = starting_index;
	}
	return false;
}

void wfa::next(wavefront_t& wavefront, int32_t s, int32_t x, int32_t o, int32_t e) {
	int32_t m_high_sx = 0;
	int32_t m_low_sx = 0;
	int32_t m_high_soe = 0;
	int32_t m_low_soe = 0;
	int32_t i_high_se = 0;
	int32_t i_low_se = 0;
	int32_t d_high_se = 0;
	int32_t d_low_se = 0;
	if (s - x >= 0) {
		m_high_sx = wavefront.wave_size(s - x, false);
		m_low_sx = wavefront.wave_size(s - x, true);
	}
	if (s - o - e >= 0) {
		m_high_soe = wavefront.wave_size(s - o - e, false);
		m_low_soe = wavefront.wave_size(s - o - e, true);
	}
	if (s - e >= 0) {
		i_high_se = wavefront.wave_size(s - e, false);
		i_low_se = wavefront.wave_size(s - e, true);
		d_high_se = wavefront.wave_size(s - e, false);
		d_low_se = wavefront.wave_size(s - e, true);
	}
	
	int32_t high = std::max({m_high_sx, m_high_soe, i_high_se, d_high_se}) + 1;
	int32_t low = std::min({m_low_sx, m_low_soe, i_low_se, d_low_se}) - 1;

	auto& low_hi_cur = wavefront.low_hi.emplace_back();
	low_hi_cur[0] = low;
	low_hi_cur[1] = high;

	auto& wave_cur = wavefront.data.emplace_back();
	wavefront.score_to_index.emplace(s, static_cast<int32_t>(wavefront.data.size()) - 1);
	for (auto& vec : wave_cur) {
		vec.resize(high - low + 1);
	}
	for (int32_t k = low; k < high + 1; ++k) {
		wave_cur[ins][k - low] = std::max({wavefront.lookup(s - o - e, match, k - 1 - low), wavefront.lookup(s - e, ins, k - 1 - low)}) + 1;
		wave_cur[del][k - low] = std::max({ wavefront.lookup(s - o - e, match, k + 1 - low), wavefront.lookup(s - e, del, k + 1 - low) });
		wave_cur[match][k - low] = std::max({ wave_cur[ins][k - low], wave_cur[del][k - low], wavefront.lookup(s - x, match, k - low) + 1});
	}
}

int32_t wfa::wavefront_simd(std::string_view a, std::string_view b, int32_t x, int32_t o, int32_t e) {
	wavefront_t wavefront;
	wavefront.data.emplace_back(std::array<std::vector<int32_t>, 3>{std::vector<int32_t>{0}, std::vector<int32_t>{0}, std::vector<int32_t>{0}});
	
	wavefront.low_hi.emplace_back(std::array{0,0});
	wavefront.score_to_index.emplace(0, 0);
	bool matched = false;

	int32_t score = 0;
	while (not matched) {
		matched = extend(wavefront, a, b, score);
		if (matched) {
			break;
		}
		if (score == 0) {
			score += std::min(x, o + e);
		}

		else {
			score = score + 1;
			while (not (
				wavefront.score_to_index.contains(score - x) 
				or wavefront.score_to_index.contains(score - e - o) 
				or (wavefront.score_to_index.contains(score - e) and (score - e >= o))
			)) {
				++score;
			}
		}
		next(wavefront, score, x, o, e);
	}
	wavefront.print();
	return score;
}