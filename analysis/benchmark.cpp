#include <iostream>
#include <string>
#include <vector>
#include <filesystem>
#include <cstdlib> // for std::system

int main() {
    if (argc < 6) {
        std::cerr << "Usage: " << argv[0] << " <seq1> <seq2> <x> <o> <e>\n";
        return 1;
    }

    std::string seq1 = argv[1];
    std::string seq2 = argv[2];
    int x = std::atoi(argv[3]);
    int o = std::atoi(argv[4]);
    int e = std::atoi(argv[5]);


}
