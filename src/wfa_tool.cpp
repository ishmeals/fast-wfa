#include "include/data_gen.hpp"
#include "include/kokkos_simd.hpp"

int main() {
    /*auto sequences = wfa::modify_sequences(100, 1000, 0.25);
    for (const auto& pair : sequences) {
        int32_t dist = wfa::wavefront_simd(pair.first, pair.second, 4, 6, 2);

    }*/

    std::string a = "tert";
    std::string b = "test";
    wfa::wavefront_simd(a, b, 4, 6, 2);
}