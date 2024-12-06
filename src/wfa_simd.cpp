#include "include/wfa_simd.hpp"

#include "fmt/format.h"
#include "fmt/ranges.h"



bool wfa::extend_simd(wavefront_t& wavefront, std::string_view a, std::string_view b, int32_t score) {
	std::vector<int32_t>& matchfront_back = wavefront.data.back()[2];
	int32_t k_low = wavefront.wave_size(score, true);
	int32_t k_high = wavefront.wave_size(score, false);
	for (int32_t k = k_low; k < k_high + 1; ++k) {
		int32_t starting_index = matchfront_back[k - k_low];
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
			//fmt::println("{}", starting_index);
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
				if (!Kokkos::Experimental::all_of(mask)) {
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
				}
				else {
					starting_index += simd_size;
					v += simd_size;
					h += simd_size;
				}
			}
		}
		
		matchfront_back[k - k_low] = starting_index;
	}
	return false;
}

void wfa::next_simd(wavefront_t& wavefront, int32_t s, int32_t x, int32_t o, int32_t e) {
	int32_t m_high_sx = wavefront.wave_size(s - x, false);
	int32_t m_low_sx = wavefront.wave_size(s - x, true);
	int32_t m_high_soe = wavefront.wave_size(s - o - e, false);
	int32_t m_low_soe = wavefront.wave_size(s - o - e, true);
	int32_t i_high_se = wavefront.wave_size(s - e, false);
	int32_t i_low_se = wavefront.wave_size(s - e, true);
	int32_t d_high_se = wavefront.wave_size(s - e, false);
	int32_t d_low_se = wavefront.wave_size(s - e, true);
	
	
	int32_t high = std::max({m_high_sx, m_high_soe, i_high_se, d_high_se}) + 1;
	int32_t low = std::min({m_low_sx, m_low_soe, i_low_se, d_low_se}) - 1;
	if (s < o + e) {
		high = 0;
		low = 0;
	}
	auto& low_hi_cur = wavefront.low_hi.emplace_back();
	low_hi_cur[0] = low;
	low_hi_cur[1] = high;

	wavefront.valid_scores.emplace(s);
	auto& wave_cur = wavefront.data.emplace_back(std::array<std::vector<int32_t>, 3>{
		std::vector<int32_t>(high - low + 1),
		std::vector<int32_t>(high - low + 1),
		std::vector<int32_t>(high - low + 1)
	});

	constexpr int32_t simd_size = static_cast<int32_t>(simd_type::size());

	int32_t k = low;
	while (high - k + 1 > simd_size) {
		auto M_soe_down = [&wavefront, s, o, e, k] (size_t i) {return wavefront.lookup(s - o - e, match, k + static_cast<int32_t>(i) - 1);};
		auto I_se_down = [&wavefront, s, e, k](size_t i) {return wavefront.lookup(s - e, ins, k + static_cast<int32_t>(i) - 1);};
		simd_type M_soe_down_vec(M_soe_down);
		simd_type I_se_down_vec(I_se_down);
		simd_type I = Kokkos::max(M_soe_down_vec, I_se_down_vec);
		Kokkos::Experimental::where(I != -1, I) = I + 1;

		auto M_soe_up = [&wavefront, s, o, e, k](size_t i) {return wavefront.lookup(s - o - e, match, k + static_cast<int32_t>(i) + 1); };
		auto D_se_up = [&wavefront, s, e, k](size_t i) {return wavefront.lookup(s - e, del, k + static_cast<int32_t>(i) + 1); };

		simd_type M_soe_up_vec(M_soe_up);
		simd_type D_se_up_vec(D_se_up);
		simd_type D = Kokkos::max(M_soe_up_vec, D_se_up_vec);

		auto M_sx = [&wavefront, s, x, k] (size_t i) {return wavefront.lookup(s - x, match, k + static_cast<int32_t>(i));};
		simd_type M_sx_vec(M_sx);
		Kokkos::Experimental::where(M_sx_vec != -1, M_sx_vec) = M_sx_vec + 1;
		simd_type M = Kokkos::max(I, D);
		M = Kokkos::max(M, M_sx_vec);

		I.copy_to(wave_cur[ins].data() + (k - low), tag_type());
		D.copy_to(wave_cur[del].data() + (k - low), tag_type());
		M.copy_to(wave_cur[match].data() + (k - low), tag_type());

		k += simd_size;
		//fmt::println("{}", simd_size);
	}
	//fmt::println("skipped to k = {}", k);
	for (/*int32_t k = low*/; k < high + 1; ++k) {
		wave_cur[ins][k - low] = std::max({wavefront.lookup(s - o - e, match, k - 1), wavefront.lookup(s - e, ins, k - 1)});
		if (wave_cur[ins][k - low] != -1) {
			++wave_cur[ins][k - low];
		}
		wave_cur[del][k - low] = std::max({ wavefront.lookup(s - o - e, match, k + 1), wavefront.lookup(s - e, del, k + 1) });
		int32_t match_sxk = wavefront.lookup(s - x, match, k);
		if (match_sxk < 0) {
			wave_cur[match][k - low] = std::max({ wave_cur[ins][k - low], wave_cur[del][k - low], -1 });
		}
		else {
			wave_cur[match][k - low] = std::max({ wave_cur[ins][k - low], wave_cur[del][k - low], match_sxk + 1});
		}
	}

}

int32_t wfa::wavefront_simd(std::string_view a, std::string_view b, int32_t x, int32_t o, int32_t e) {
	/*std::vector<int32_t> a_int;
	std::vector<int32_t> b_int;
	for (const char c : a) {
		a_int.emplace_back(static_cast<int32_t>(c));
	}
	for (const char c : b) {
		b_int.emplace_back(static_cast<int32_t>(c));
	}*/

	wavefront_t wavefront;
	wavefront.data.emplace_back(std::array<std::vector<int32_t>, 3>{std::vector<int32_t>{-1}, std::vector<int32_t>{-1}, std::vector<int32_t>{0}});
	
	wavefront.low_hi.emplace_back(std::array{0,0});
	wavefront.valid_scores.emplace(0);
	bool matched = false;

	int32_t score = 0;

	int32_t final_k = static_cast<int32_t>(b.size()) - static_cast<int32_t>(a.size());
	int32_t final_offset = static_cast<int32_t>(b.size());

	while (not matched) {
		extend_simd(wavefront, a, b, score);
		std::vector<int32_t>& matchfront_back = wavefront.data.back()[2];
		int32_t k_low = wavefront.wave_size(score, true);
		if (matchfront_back[final_k - k_low] >= final_offset) {
			break;
		}
		else {
			score = score + 1;
			while (not (
				wavefront.valid_scores.contains(score - x)
				or wavefront.valid_scores.contains(score - e - o)
				or (wavefront.valid_scores.contains(score - e) and (score - e >= o))
			)) {
				wavefront.data.emplace_back(std::array<std::vector<int32_t>, 3>{std::vector<int32_t>{-1}, std::vector<int32_t>{-1}, std::vector<int32_t>{-1}});
				wavefront.low_hi.emplace_back(std::array{ 0,0 });
				++score;
			}
		}
		next_simd(wavefront, score, x, o, e);
		//wavefront.print();
		//fmt::println("---");
	}
	//wavefront.print();
	return -1*score;
}