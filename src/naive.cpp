#include "include/naive.hpp"
#include <algorithm>
#include <array>
#include <string_view>
#include <vector>
#include <limits>

 int32_t wfa::naive(std::string_view a, std::string_view b,  int32_t x,  int32_t o,  int32_t e) {
    constexpr int32_t INF = static_cast<int32_t>(-1e9);
    size_t n = a.size(), m = b.size();
    std::vector<std::vector<std::array<int32_t, 3>>> dp(n + 1, std::vector<std::array<int, 3>>(m + 1, {INF, INF, INF})); // M, I, D

    dp[0][0][0] = 0; // Initial state at (0, 0)

    // First column (gaps in B)
    for (size_t i = 1; i <= n; ++i) {
        dp[i][0][0] = -o - e * static_cast<int32_t>(i);  // Deletion costs
        dp[i][0][2] = dp[i][0][0];
    }
    // First row (gaps in A)
    for (size_t j = 1; j <= m; ++j) {
        dp[0][j][0] = -o - e * static_cast<int32_t>(j);  // Insertion costs
        dp[0][j][1] = dp[0][j][0];
    }

    // DP computation
    for (size_t i = 1; i <= n; i++) {
        for (size_t j = 1; j <= m; j++) {
            int32_t cost = (a[i - 1] != b[j - 1]) ? -x : 0;
            dp[i][j][0] = std::max({dp[i - 1][j - 1][0], dp[i - 1][j - 1][1], dp[i - 1][j - 1][2]}) + cost; // Match/Mismatch state
            dp[i][j][1] = std::max(dp[i][j - 1][0] - o - e, dp[i][j - 1][1] - e); // Insertion state (gap in A)
            dp[i][j][2] = std::max(dp[i - 1][j][0] - o - e, dp[i - 1][j][2] - e); // Deletion state (gap in B)
        }
    }

    // Final result
    return std::max({dp[n][m][0], dp[n][m][1], dp[n][m][2]});
}

 int32_t wfa::wavefront_dp(std::string_view a, std::string_view b,  int32_t x,  int32_t o,  int32_t e) {
     constexpr int32_t INF = static_cast<int32_t>(-1e9);  // Use large negative value to avoid overflow
    size_t n = a.size(), m = b.size();
    std::vector<std::vector<std::array<int, 3>>> dp(
        n + 1, std::vector<std::array<int, 3>>(m + 1, {INF, INF, INF})); // M, I, D
    
    dp[0][0] = {0, INF, INF};  // Initialize starting point

    // First row (gaps in A)
    for (size_t j = 1; j <= m; ++j) {
        dp[0][j][1] = -o - e * static_cast<int32_t>(j);  // Insertion costs
        dp[0][j][0] = dp[0][j][1];
    }

    // First column (gaps in B)
    for (size_t i = 1; i <= n; ++i) {
        dp[i][0][2] = -o - e * static_cast<int32_t>(i);  // Deletion costs
        dp[i][0][0] = dp[i][0][2];
    }

    // Queue to manage wavefront processing at each cost level
    std::vector<std::vector<std::pair<size_t, size_t>>> q(1);
    q[0].emplace_back(0, 0);

    // Add initial positions from the first row and column to the queue
    for (size_t i = 1; i <= n; ++i) {
         int32_t cost = dp[i][0][0];
        if (q.size() <= static_cast<size_t>(-cost)) q.resize(-cost + 1);
        q[-cost].emplace_back(i, 0);
    }
    for (size_t j = 1; j <= m; ++j) {
         int32_t cost = dp[0][j][0];
        if (q.size() <= static_cast<size_t>(-cost)) q.resize(-cost + 1);
        q[-cost].emplace_back(0, j);
    }

    // Process cells in wavefront order
    for (size_t s = 0; s < q.size(); ++s) {
        for (size_t k = 0; k < q[s].size(); ++k) {
            auto [i, j] = q[s][k];

             int32_t max_dp_cost = std::max({dp[i][j][0], dp[i][j][1], dp[i][j][2]});
            if (s > static_cast<size_t>(-max_dp_cost)) continue;

            if (i == n && j == m) {
                return max_dp_cost;
            }

            // Match/Mismatch case
            if (i < n && j < m) {
                 int32_t match_cost = (a[i] != b[j]) ? -x : 0;
                 int32_t prev_max = std::max({dp[i][j][0], dp[i][j][1], dp[i][j][2]});
                 int32_t cost = prev_max + match_cost;
                if (cost > dp[i + 1][j + 1][0]) {
                    dp[i + 1][j + 1][0] = cost;
                    if (q.size() <= static_cast<size_t>(-cost)) q.resize(-cost + 1);
                    q[-cost].emplace_back(i + 1, j + 1);
                }
            }

            // Insertion (gap in A)
            if (j < m) {
                 int32_t insert_cost_new = dp[i][j][0] - o - e;
                 int32_t insert_cost_extend = dp[i][j][1] - e;
                 int32_t cost = std::max(insert_cost_new, insert_cost_extend);
                if (cost > dp[i][j + 1][1]) {
                    dp[i][j + 1][1] = cost;
                    if (q.size() <= static_cast<size_t>(-cost)) q.resize(-cost + 1);
                    q[-cost].emplace_back(i, j + 1);
                }
            }

            // Deletion (gap in B)
            if (i < n) {
                 int32_t delete_cost_new = dp[i][j][0] - o - e;
                 int32_t delete_cost_extend = dp[i][j][2] - e;
                 int32_t cost = std::max(delete_cost_new, delete_cost_extend);
                if (cost > dp[i + 1][j][2]) {
                    dp[i + 1][j][2] = cost;
                    if (q.size() <= static_cast<size_t>(-cost)) q.resize(-cost + 1);
                    q[-cost].emplace_back(i + 1, j);
                }
            }
        }
    }

    return -INF; // Return the minimum cost
}
