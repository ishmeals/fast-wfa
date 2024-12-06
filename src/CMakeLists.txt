set(WFA_SRC_FILES ${WFA_SRC_FILES}
	"src/naive.cpp"
	"src/wfa_simd.cpp"
	"src/wfa.cpp"
	"src/data_gen.cpp"
	"src/parallel.cpp"
	PARENT_SCOPE
)

set(WFA_HEADER_FILES ${WFA_HEADER_FILES}
	"include/naive.hpp"
	"include/wfa_simd.hpp"
	"include/wfa.hpp"
	"include/data_gen.hpp"
	"include/parallel.hpp"
	PARENT_SCOPE
)