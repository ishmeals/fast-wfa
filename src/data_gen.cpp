#include "include/data_gen.hpp"
#include <random>
#include <array>

constexpr std::array<char, 4> nt = { 'A', 'T', 'G', 'C' };

std::vector<std::string> wfa::generate_sequences(int32_t length, int32_t count) {
	std::vector<std::string> strings;

	std::random_device r;
	std::default_random_engine random(r());
	std::uniform_int_distribution<int32_t> int_rand(0, 3);

	for (int32_t i = 0; i < count; ++i) {
		std::string str;
		for (int32_t j = 0; j < length; ++j) {
			str.push_back(nt[int_rand(random)]);
		}
		strings.emplace_back(std::move(str));
	}

	return strings;
}

std::vector<std::pair<std::string, std::string>> wfa::modify_sequences(int32_t length, int32_t count, double error_rate) {
	auto strings = generate_sequences(length, count);

	std::random_device r;
	std::default_random_engine random(r());
	std::uniform_int_distribution<int32_t> int_rand(0, 3);
	std::uniform_real_distribution<double> real_rand(0, 1);

	std::vector<std::pair<std::string, std::string>> sequences;
	for (auto& str : strings) {
		std::string modified = str;
		for (char& c : modified) {
			if (real_rand(random) < error_rate) {
				c = nt[int_rand(random)];
			}
		}
		sequences.emplace_back(std::pair{str, modified});
	}

	return sequences;
}