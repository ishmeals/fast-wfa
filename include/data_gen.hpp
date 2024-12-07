#pragma once

#include <vector>
#include <string>


namespace wfa {
	//generates sequences of nucleobases of provided length and count
	std::vector<std::string> generate_sequences(int32_t length, int32_t count);

	//generates pairs sequences of nucleobases of provided length and count, where the second element in the pair the first sequence modified by the provided error rate.
	std::vector<std::pair<std::string, std::string>> modify_sequences(int32_t length, int32_t count, double error_rate);

}