#include "include/data_gen.hpp"
#include "include/kokkos_simd.hpp"
#include "include/naive.hpp"
#include "fmt/format.h"


int main() {
    /*auto sequences = wfa::modify_sequences(100, 1000, 0.25);
    for (const auto& pair : sequences) {
        int32_t dist = wfa::wavefront_simd(pair.first, pair.second, 4, 6, 2);

    }*/

    std::string a = "AGGATGCTCG";
    std::string b = "ACCATACTCG";
    auto score = wfa::wavefront_simd(a, b, 4, 6, 2);
    auto score2 = wfa::naive(a, b, 4, 6, 2);
    auto score3 = wfa::wavefront(a, b, 4, 6, 2);
    fmt::println("{} {} {}", score, score2, score3);
}