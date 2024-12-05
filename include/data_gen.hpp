#pragma once

#include <vector>
#include <string>


namespace wfa {
	
	std::vector<std::string> generate_sequences(int32_t length, int32_t count);

	std::vector<std::pair<std::string, std::string>> modify_sequences(int32_t length, int32_t count, double error_rate);

}