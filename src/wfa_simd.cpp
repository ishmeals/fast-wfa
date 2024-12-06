#include "include/wfa_simd.hpp"

#include "fmt/format.h"
#include "fmt/ranges.h"



bool wfa::extend_simd(wavefront_t& wavefront, std::string_view a, std::string_view b) {
	wavefront_entry_t& entry = wavefront.views.back();
	std::span<int32_t> matchfront_back(entry.data.data() + match * entry.number_per_col, entry.number_per_col);
	for (int32_t k = entry.low; k < entry.high + 1; ++k) {
		int32_t starting_index = matchfront_back[k - entry.low];
		if (starting_index == -1) {
			continue;
		}
		int32_t v = starting_index - k;
		int32_t h = starting_index;
		
		int32_t a_size = static_cast<int32_t>(a.size());
		int32_t b_size = static_cast<int32_t>(b.size());

		bool mismatch = false;

		constexpr int32_t char_per_int = (32 / 8);
		constexpr int32_t simd_size_true = static_cast<int32_t>(simd_type::size());
		constexpr int32_t simd_size = static_cast<int32_t>(simd_type::size()) * char_per_int;

		while (not mismatch) {
			int32_t v_rem = a_size - v - 1;
			int32_t h_rem = b_size - h - 1;
			if (v_rem < simd_size or h_rem < simd_size) {
				if (v >= a_size or h >= b_size) {
					break;
				}
				const int32_t v_c = a[v];
				const int32_t h_c = b[h];
				if (v_c == h_c) {
					++starting_index;
					++v;
					++h;
				}
				else {
					mismatch = true;
				}
			}
			else {
				simd_type v_vec;
				v_vec.copy_from(reinterpret_cast<const int32_t*>(a.data() + v), tag_type());
				simd_type h_vec;
				h_vec.copy_from(reinterpret_cast<const int32_t*>(b.data() + h), tag_type());
				auto mask = v_vec == h_vec;
				//if (!Kokkos::Experimental::all_of(mask)) {
					for (int32_t i = 0; i < simd_size / char_per_int; ++i) {
						if (!mask[i]) {
							mismatch = true;
							//starting_index += i * char_per_int;
							int32_t base_offset = i * char_per_int;
							for (int32_t i = 0; i < 4; ++i) {
								if (a[v + base_offset + i] != b[h + base_offset + i]) {
									starting_index += base_offset + i;
									break;
								}
							}
							break;
						}
					}
				//}
				//else {
				if (!mismatch) {
					starting_index += simd_size;
					v += simd_size;
					h += simd_size;
				}
			}
		}
		
		matchfront_back[k - entry.low] = starting_index;
	}
	return false;
}

//return soe.lookup(match, k + static_cast<int32_t>(i) - 1);
wfa::simd_type wfa::simd_lookup(wavefront_entry_t& t, int32_t column, int32_t k, int32_t scaling) {
	constexpr int32_t simd_size = static_cast<int32_t>(simd_type::size());
	constexpr std::array<int32_t, 16> offset_arr = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15 };
	simd_type offsets;
	offsets.copy_from(offset_arr.data(), tag_type());
	simd_type const_k(k + scaling);
	auto k_range = offsets + const_k;
	auto row = k_range - t.low;
	auto bounds_mask = row >= 0 && row < t.number_per_col;
	simd_type values;

	Kokkos::Experimental::where(bounds_mask, values).copy_from(t.data.data() + (column * t.number_per_col + row[0]), tag_type());

	Kokkos::Experimental::where(not bounds_mask, values) = -1;
	
	return values;
}

void wfa::next_simd(wavefront_t& wavefront, int32_t s, int32_t x, int32_t o, int32_t e) {
	int32_t m_high_sx = wavefront.wave_size_high(s - x);
	int32_t m_low_sx = wavefront.wave_size_low(s - x);
	int32_t m_high_soe = wavefront.wave_size_high(s - o - e);
	int32_t m_low_soe = wavefront.wave_size_low(s - o - e);
	int32_t i_high_se = wavefront.wave_size_high(s - e);
	int32_t i_low_se = wavefront.wave_size_low(s - e);
	int32_t d_high_se = wavefront.wave_size_high(s - e);
	int32_t d_low_se = wavefront.wave_size_low(s - e);
	
	
	int32_t high = std::max({m_high_sx, m_high_soe, i_high_se, d_high_se}) + 1;
	int32_t low = std::min({m_low_sx, m_low_soe, i_low_se, d_low_se}) - 1;
	if (s < o + e) {
		high = 0;
		low = 0;
	}

	auto& wave_cur = wavefront.insert(low, high);

	constexpr int32_t simd_size = static_cast<int32_t>(simd_type::size());

	int32_t k = low;
	while (high - k + 1 > simd_size) {
		simd_type M_soe_down_vec;
		simd_type M_soe_up_vec;
		
		if (s - o - e >= 0) {
			int32_t soe_index = wavefront.mapping[s - o - e];
			if (soe_index != -1) {
				auto& soe = wavefront.views[soe_index];
				M_soe_down_vec = simd_lookup(soe, match, k, -1);
				M_soe_up_vec = simd_lookup(soe, match, k, +1);
			}
			else {
				M_soe_down_vec = simd_type(-1);
				M_soe_up_vec = simd_type(-1);
			}
		}
		else {
			M_soe_down_vec = simd_type(-1);
			M_soe_up_vec = simd_type(-1);
		}

		simd_type I_se_down_vec;
		simd_type D_se_up_vec;
		if (s - e >= 0) {
			int32_t se_index = wavefront.mapping[s - e];
			if (se_index != -1) {
				auto& se = wavefront.views[se_index];
				I_se_down_vec = simd_lookup(se, ins, k, -1);
				D_se_up_vec = simd_lookup(se, del, k, +1);
			}
			else {
				I_se_down_vec = simd_type(-1);
				D_se_up_vec = simd_type(-1);
			}
		}
		else {
			I_se_down_vec = simd_type(-1);
			D_se_up_vec = simd_type(-1);
		}

		simd_type M_sx_vec;
		if (s - x >= 0) {
			int32_t sx_index = wavefront.mapping[s - x];
			if (sx_index != -1) {
				auto& sx = wavefront.views[sx_index];
				M_sx_vec = simd_lookup(sx, match, k, 0);
			}
			else {
				M_sx_vec = simd_type(-1);
			}
		}
		else {
			M_sx_vec = simd_type(-1);
		}

		simd_type I = Kokkos::max(M_soe_down_vec, I_se_down_vec);
		Kokkos::Experimental::where(I != -1, I) = I + 1;

		simd_type D = Kokkos::max(M_soe_up_vec, D_se_up_vec);
		
		Kokkos::Experimental::where(M_sx_vec != -1, M_sx_vec) = M_sx_vec + 1;
		simd_type M = Kokkos::max(I, D);
		M = Kokkos::max(M, M_sx_vec);

		I.copy_to(wave_cur.start_ptr(ins) + (k - low), tag_type());
		D.copy_to(wave_cur.start_ptr(del) + (k - low), tag_type());
		M.copy_to(wave_cur.start_ptr(match) + (k - low), tag_type());

		k += simd_size;
	}

	for (; k < high + 1; ++k) {
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

int32_t wfa::wavefront_simd(std::string_view a, std::string_view b, int32_t x, int32_t o, int32_t e, wavefront_arena_t& arena) {
	wavefront_t wavefront(arena);
	auto& first = wavefront.insert(0, 0);
	first.data[0] = -1;
	first.data[1] = -1;
	first.data[2] = 0;

	bool matched = false;

	int32_t score = 0;

	int32_t final_k = static_cast<int32_t>(b.size()) - static_cast<int32_t>(a.size());
	int32_t final_offset = static_cast<int32_t>(b.size());

	while (not matched) {
		extend_simd(wavefront, a, b);
		auto& matchfront_back = wavefront.views.back();
		if (matchfront_back.lookup(match, final_k) >= final_offset) {
			break;
		}
		else {
			score = score + 1;
			while (not (
				wavefront.valid_score(score - x)
				or wavefront.valid_score(score - e - o)
				or (wavefront.valid_score(score - e) and (score - e >= o))
			)) {
				wavefront.insert();
				++score;
			}
		}
		next_simd(wavefront, score, x, o, e);
	}
	arena.current_index = 0;
	return -1*score;
}