#include "naive.hpp"
#include <algorithm>
#include <array>
#include <chrono>
#include <iostream>
#include <string>
#include <string_view>
#include <vector>

int naive(std::string_view a, std::string_view b, int x, int o, int e) {
	std::vector dp(a.size() + 1, std::vector(b.size() + 1, std::array<int,3>())); // M, I, D

	dp[0][0] = {0, o, o};
	for (int i = 1; i <= int(a.size()); i++) {
		dp[0][i] = dp[i][0] = {100'000'000, o + i * e, o + i * e};
	}
	for (size_t i = 1; i <= a.size(); i++) {
		for (size_t j = 1; j <= b.size(); j++) {
			dp[i][j][1] = std::min(dp[i][j-1][0] + o + e, dp[i][j-1][1] + e);
			dp[i][j][2] = std::min(dp[i-1][j][0] + o + e, dp[i-1][j][2] + e);
			dp[i][j][0] = std::min({dp[i-1][j-1][0] + (a[i-1] != b[j-1]) * x, dp[i][j][1], dp[i][j][2]});
		}
	}
	return dp.back().back()[0];
}

int wavefront(std::string_view a, std::string_view b, int x, int o, int e) {
	std::vector dp(a.size() + 1, std::vector(b.size() + 1, std::array<int,3>({100'000'000, 100'000'000, 100'000'000}))); // M, I, D

	dp[0][0] = {0, o, o};
	std::vector<std::vector<std::pair<size_t,size_t>>> q(1, {{0,0}});
	for (int s = 0; ; s++) {
		for (size_t k = 0; k < q[s].size(); k++) {
			auto [i, j] = q[s][k];
			if (i == a.size() && j == b.size()) return s;

			if (i < a.size() && j < b.size()) {
				int cost = dp[i][j][0] + (a[i] != b[j]) * x;
				if (cost < dp[i+1][j+1][0]) {
					dp[i+1][j+1][0] = cost;
					while (q.size() <= cost) q.push_back({});
					q[cost].push_back({i+1, j+1});
				}
			}
			if (i < a.size()) {
				int cost = std::min(dp[i][j][0] + o + e, dp[i][j][2] + e);
				if (cost < dp[i+1][j][2]) {
					dp[i+1][j][2] = cost;
					dp[i+1][j][0] = std::min(dp[i+1][j][0], cost);
					while (q.size() <= cost) q.push_back({});
					q[cost].push_back({i+1, j});
				}
			} 
			if (j < b.size()) {
				int cost = std::min(dp[i][j][0] + o + e, dp[i][j][1] + e);
				if (cost < dp[i][j+1][1]) {
					dp[i][j+1][1] = cost;
					dp[i][j+1][0] = std::min(dp[i][j+1][0], cost);
					while (q.size() <= cost) q.push_back({});
					q[cost].push_back({i, j+1});
				}
			}
		}
	}
}


int main() {
	/*std::string a("GATACA"), b("GAGATA");*/

	const int N = 2000;
	std::string a (N, ' '), b(N, ' ');
	for (int i = 0; i < N; i++) {
		a[i] = b[i] = "ACGT"[rand() % 4];
		if (rand() % 20) continue;
		while (a[i] == b[i]) a[i] =  "ACGT"[rand() % 4];
	}

	auto t = std::chrono::system_clock::now();
	std::cout << naive(a,b, 4, 6, 2) << std::endl;
	std::cout << std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::now() - t).count() << std::endl;

	t = std::chrono::system_clock::now();
	std::cout << wavefront(a,b, 4, 6, 2) << std::endl;
	std::cout << std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::now() - t).count() << std::endl;

}
