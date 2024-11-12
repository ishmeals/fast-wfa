#include "naive.hpp"
#include <algorithm>
#include <array>
#include <chrono>
#include <iostream>
#include <string>
#include <string_view>
#include <vector>

int naive(std::string_view a, std::string_view b, int x, int o, int e) {
    const int INF = 1e9;
    std::vector dp(a.size() + 1, std::vector(b.size() + 1, std::array<int, 3>{INF, INF, INF})); // M, I, D

    dp[0][0] = {0, INF, INF}; // Initial state at (0, 0)

    // First column (gaps in B)
    for (size_t i = 1; i <= a.size(); i++) {
        dp[i][0][0] = dp[i][0][1] = INF;
        dp[i][0][2] = std::min(dp[i - 1][0][0] + o, dp[i - 1][0][2] + e);
    }
    // First row (gaps in A)
    for (size_t j = 1; j <= b.size(); j++) {
        dp[0][j][0] = dp[0][j][2] = INF;
        dp[0][j][1] = std::min(dp[0][j - 1][0] + o, dp[0][j - 1][1] + e);
    }

    // DP computation
    for (size_t i = 1; i <= a.size(); i++) {
        for (size_t j = 1; j <= b.size(); j++) {
            int cost = (a[i - 1] != b[j - 1]) * x;

            // Match/Mismatch state
            dp[i][j][0] = std::min({
                dp[i - 1][j - 1][0] + cost,
                dp[i - 1][j - 1][1] + cost,
                dp[i - 1][j - 1][2] + cost
            });

            dp[i][j][1] = std::min(dp[i][j - 1][0] + o, dp[i][j - 1][1] + e); // Insertion state (gap in A)
            dp[i][j][2] = std::min(dp[i - 1][j][0] + o, dp[i - 1][j][2] + e); // Deletion state (gap in B)
        }
    }

    // Final result
    return std::min({dp[a.size()][b.size()][0], dp[a.size()][b.size()][1], dp[a.size()][b.size()][2]});
}


int wavefront(std::string_view a, std::string_view b, int x, int o, int e) {
    const int INF = 1e9;
    std::vector dp(a.size() + 1, std::vector(b.size() + 1, std::array<int, 3>{INF, INF, INF})); // M, I, D

    dp[0][0] = {0, INF, INF}; // Initialize starting point

    // First row (gaps in A)
    for (size_t j = 1; j <= b.size(); ++j) {
        dp[0][j][1] = o + e * (j - 1);
        dp[0][j][0] = dp[0][j][1];
    }
    // First column (gaps in B)
    for (size_t i = 1; i <= a.size(); ++i) {
        dp[i][0][2] = o + e * (i - 1);
        dp[i][0][0] = dp[i][0][2];
    }

    // Queue to manage wavefront processing at each cost level
    std::vector<std::vector<std::pair<size_t, size_t>>> q(1, {{0, 0}});

    // Process cells in wavefront order
    for (int s = 0; s < q.size(); ++s) {
        for (size_t k = 0; k < q[s].size(); ++k) {
            auto [i, j] = q[s][k];
            if (i == a.size() && j == b.size()) return s;

            // Match/Mismatch case
            if (i < a.size() && j < b.size()) {
                int cost = dp[i][j][0] + (a[i] != b[j]) * x;
                if (cost < dp[i + 1][j + 1][0]) {
                    dp[i + 1][j + 1][0] = cost;
                    while (q.size() <= cost) q.push_back({});
                    q[cost].emplace_back(i + 1, j + 1);
                }
            }

            // Insertion (gap in A)
            if (j < b.size()) {
                int cost = std::min(dp[i][j][0] + o, dp[i][j][1] + e);
                if (cost < dp[i][j + 1][1]) {
                    dp[i][j + 1][1] = cost;
                    dp[i][j + 1][0] = std::min(dp[i][j + 1][0], cost);
                    while (q.size() <= cost) q.push_back({});
                    q[cost].emplace_back(i, j + 1);
                }
            }

            // Deletion (gap in B)
            if (i < a.size()) {
                int cost = std::min(dp[i][j][0] + o, dp[i][j][2] + e);
                if (cost < dp[i + 1][j][2]) {
                    dp[i + 1][j][2] = cost;
                    dp[i + 1][j][0] = std::min(dp[i + 1][j][0], cost);
                    while (q.size() <= cost) q.push_back({});
                    q[cost].emplace_back(i + 1, j);
                }
            }
        }
    }

    return dp[a.size()][b.size()][0]; // Return the minimum cost at the bottom-right corner
}

/**
int main() {
	/*std::string a("GATACA"), b("GAGATA");*/ /*

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
*/