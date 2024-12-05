set(WFA_SRC_FILES ${WFA_SRC_FILES}
	"src/naive.cpp"
	"src/kokkos_simd.cpp"
	"src/data_gen.cpp"
	PARENT_SCOPE
)

set(WFA_HEADER_FILES ${WFA_HEADER_FILES}
	"include/naive.hpp"
	"include/kokkos_simd.hpp"
	"include/data_gen.hpp"
	PARENT_SCOPE
)