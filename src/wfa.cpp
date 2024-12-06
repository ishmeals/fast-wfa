#include "include/wfa.hpp"

#include "fmt/format.h"
#include "fmt/ranges.h"
#include <span>

int32_t wfa::wavefront_entry_t::lookup(int32_t column, int32_t k) {
	int32_t row = k - low;
	if (row < 0) {
		return -1;
	}
	if (row > high - low) {
		return -1;
	}

	return data[column * number_per_col + row];
}

void wfa::wavefront_entry_t::set(int32_t column, int32_t k, int32_t val) {
	int32_t row = k - low;
	data[column * number_per_col + row] = val;
}

int32_t& wfa::wavefront_entry_t::no_bound(int32_t column, int32_t row) {
	return data[column * number_per_col + row];
}

int32_t* wfa::wavefront_entry_t::start_ptr(int32_t column) {
	return data.data() + column * number_per_col;
}

wfa::wavefront_entry_t::wavefront_entry_t(int32_t low, int32_t high) : low(low), high(high), number_per_col(high - low + 1) {
	data.resize(3 * number_per_col);
}

int32_t wfa::wavefront_t::lookup(int32_t score, int32_t column, int32_t k) {
	if (score < 0) {
		return -1;
	}
	/*auto iter = score_to_index.find(score);
	if (iter == score_to_index.end()) {
		return -1;
	}
	int32_t low_correction = wave_size(score, true);
	int32_t row = k - low_correction;
	if (row < 0) {
		return -1;
	}
	//if (row > static_cast<int32_t>(data[iter->second][column].size()) - 1) {
	if (row > static_cast<int32_t>(data[score][column].size()) - 1) {
		return -1;
	}
	return (data[score][column][row]);*/
	return data[score].lookup(column, k);
}

int32_t wfa::wavefront_t::wave_size_low(int32_t score) {
	/*auto iter = score_to_index.find(score);
	if (iter == score_to_index.end()) {
		return 0;
	}
	if (score < 0) {
		return -1;
	}

	auto& pair = low_hi[score];
	if (low) {
		return pair[0];
	}
	return pair[1]*/;
	if (score < 0) {
		return -1;
	}
	return data[score].low;
}

int32_t wfa::wavefront_t::wave_size_high(int32_t score) {
	/*auto iter = score_to_index.find(score);
	if (iter == score_to_index.end()) {
		return 0;
	}
	if (score < 0) {
		return -1;
	}

	auto& pair = low_hi[score];
	if (low) {
		return pair[0];
	}
	return pair[1]*/;
	if (score < 0) {
		return -1;
	}
	return data[score].high;
}




void wfa::wavefront_t::print() {
	/*for (const auto& pair : score_to_index) {
		fmt::println("Score: {}", pair.first);
		fmt::println("I: {}\nD: {}\nM: {}", data[pair.second][ins], data[pair.second][del], data[pair.second][match]);
	}*/
}

bool wfa::extend(wavefront_t& wavefront, std::string_view a, std::string_view b) {
	//std::vector<int32_t>& matchfront_back = wavefront.data.back()[2];
	//int32_t k_low = wavefront.wave_size(score, true);
	//int32_t k_high = wavefront.wave_size(score, false);
	wavefront_entry_t& entry = wavefront.data.back();
	std::span<int32_t> matchfront_back(entry.data.data() + match * entry.number_per_col, entry.number_per_col);
	for (int32_t k = entry.low; k < entry.high + 1; ++k) {
		int32_t starting_index = matchfront_back[k - entry.low];
		if (starting_index == -1) {
			continue;
		}
		int32_t v = starting_index - k;
		int32_t h = starting_index;
		bool mismatch = false;
		while (not mismatch) {
			if (v >= static_cast<int32_t>(a.size()) or h >= static_cast<int32_t>(b.size())) {
				//fmt::println("Extend out of bounds with: {} {}", v, h);
				break;
			}
			const char v_c = a[v];
			const char h_c = b[h];
			if (v_c == h_c) {
				/*if (v == static_cast<int32_t>(a.size()) - 1 and h == static_cast<int32_t>(b.size()) - 1) {
					return true;
				}*/
				++starting_index;
				++v;
				++h;
			}
			else {
				mismatch = true;
			}
		}
		matchfront_back[k - entry.low] = starting_index;
	}
	return false;
}

void wfa::next(wavefront_t& wavefront, int32_t s, int32_t x, int32_t o, int32_t e) {
	/*int32_t m_high_sx = 0;
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
	}*/

	int32_t m_high_sx = wavefront.wave_size_high(s - x);
	int32_t m_low_sx = wavefront.wave_size_low(s - x);
	int32_t m_high_soe = wavefront.wave_size_high(s - o - e);
	int32_t m_low_soe = wavefront.wave_size_low(s - o - e);
	int32_t i_high_se = wavefront.wave_size_high(s - e);
	int32_t i_low_se = wavefront.wave_size_low(s - e);
	int32_t d_high_se = wavefront.wave_size_high(s - e);
	int32_t d_low_se = wavefront.wave_size_low(s - e);


	int32_t high = std::max({ m_high_sx, m_high_soe, i_high_se, d_high_se }) + 1;
	int32_t low = std::min({ m_low_sx, m_low_soe, i_low_se, d_low_se }) - 1;
	if (s < o + e) {
		high = 0;
		low = 0;
	}
	/*auto& low_hi_cur = wavefront.low_hi.emplace_back();
	low_hi_cur[0] = low;
	low_hi_cur[1] = high;*/

	auto& wave_cur = wavefront.data.emplace_back(low, high);

	wavefront.valid_scores.emplace(s);
	/*auto& wave_cur = wavefront.data.emplace_back(std::array<std::vector<int32_t>, 3>{
		std::vector<int32_t>(high - low + 1),
			std::vector<int32_t>(high - low + 1),
			std::vector<int32_t>(high - low + 1)
	});*/

	for (int32_t k = low; k < high + 1; ++k) {
		wave_cur.no_bound(ins, k - low) = std::max({ wavefront.lookup(s - o - e, match, k - 1), wavefront.lookup(s - e, ins, k - 1) });
		if (wave_cur.no_bound(ins, k - low) != -1) {
			++wave_cur.no_bound(ins, k - low);
		}
		wave_cur.no_bound(del, k - low) = std::max({ wavefront.lookup(s - o - e, match, k + 1), wavefront.lookup(s - e, del, k + 1) });
		int32_t match_sxk = wavefront.lookup(s - x, match, k);
		if (match_sxk < 0) {
			wave_cur.no_bound(match, k - low) = std::max({ wave_cur.no_bound(ins, k - low), wave_cur.no_bound(del, k - low), -1 });
		}
		else {
			wave_cur.no_bound(match, k - low) = std::max({ wave_cur.no_bound(ins, k - low), wave_cur.no_bound(del, k - low), match_sxk + 1 });
		}

	}
}

int32_t wfa::wavefront(std::string_view a, std::string_view b, int32_t x, int32_t o, int32_t e) {
	wavefront_t wavefront;
	auto& first = wavefront.data.emplace_back(0,0);
	first.data = {-1, -1, 0};

	//wavefront.low_hi.emplace_back(std::array{ 0,0 });
	wavefront.valid_scores.emplace(0);
	bool matched = false;

	int32_t score = 0;

	int32_t final_k = static_cast<int32_t>(b.size()) - static_cast<int32_t>(a.size());
	int32_t final_offset = static_cast<int32_t>(b.size());

	while (not matched) {
		/*matched =*/ extend(wavefront, a, b);
		auto& matchfront_back = wavefront.data.back();
		//int32_t k_low = wavefront.wave_size(score, true);
		if (matchfront_back.lookup(match, final_k) >= final_offset) {
			break;
		}
		else {
			score = score + 1;
			while (not (
				wavefront.valid_scores.contains(score - x)
				or wavefront.valid_scores.contains(score - e - o)
				or (wavefront.valid_scores.contains(score - e) and (score - e >= o))
				)) {
				//wavefront.data.emplace_back(std::array<std::vector<int32_t>, 3>{std::vector<int32_t>{-1}, std::vector<int32_t>{-1}, std::vector<int32_t>{-1}});
				wavefront.data.emplace_back(0,0);
				//wavefront.low_hi.emplace_back(std::array{ 0,0 });
				++score;
			}
		}
		next(wavefront, score, x, o, e);
		//wavefront.print();
		//fmt::println("---");
	}
	//wavefront.print();
	return -1 * score;
}