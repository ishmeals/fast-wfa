cmake_minimum_required (VERSION 3.28)

include(FetchContent)

#set(VCPKG_OVERLAY_TRIPLETS ${CMAKE_SOURCE_DIR}/vcpkg_overlays/ ${CMAKE_SOURCE_DIR}/my-vcpkg-triplets)
#set(CMAKE_TOOLCHAIN_FILE ${CMAKE_SOURCE_DIR}/vcpkg/scripts/buildsystems/vcpkg.cmake)
#set(CMAKE_EXPORT_COMPILE_COMMANDS 1)
# if(CMAKE_CXX_COMPILER_ID STREQUAL "Clang")
	# set(CMAKE_LINKER_TYPE LLD)
# endif()
#Static Analysis
#find_program(CMAKE_CXX_CPPCHECK cppcheck NAMES cppcheck HINTS $ENV{PROGRAMFILES}/cppcheck)
#if(CMAKE_CXX_CPPCHECK)
#	list(APPEND CMAKE_CXX_CPPCHECK "--enable=all" "--suppressions-list=${CMAKE_SOURCE_DIR}/suppressions.txt" "--template=vs")
#endif()

project ("SynthesisSearch")


#if(CMAKE_CXX_COMPILER_ID STREQUAL "Clang")
#	set(CMAKE_LINKER_TYPE LLD)
#endif()

set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)
set(CMAKE_CXX_STANDARD 20) # Needed for Kokkos

#
# Dependencies
#

FetchContent_Declare(
	kokkos 
	GIT_REPOSITORY https://github.com/kokkos/kokkos
	GIT_TAG 15dc143e5f39949eece972a798e175c4b463d4b8
	GIT_SHALLOW TRUE
)

FetchContent_Declare(
  fmt
  GIT_REPOSITORY https://github.com/fmtlib/fmt
  GIT_TAG        e69e5f977d458f2650bb346dadf2ad30c5320281
  GIT_SHALLOW TRUE
) # 10.2.1

# Fetch WFA2-lib
FetchContent_Declare(
  WFA2lib
  GIT_REPOSITORY https://github.com/QuantumFelidae/WFA2-lib.git
  GIT_TAG        main
)


set(ENABLE_TESTING OFF)
FetchContent_MakeAvailable(fmt kokkos WFA2lib)

add_subdirectory(src)

#
# LIB_WFA
#

set(LIB_WFA ${LIB_WFA} fast_wfa_lib)
add_library(${LIB_WFA} STATIC)

target_sources(${LIB_WFA} 
	PUBLIC FILE_SET headers TYPE HEADERS FILES
	${WFA_HEADER_FILES}
)

target_sources(${LIB_WFA} PRIVATE ${WFA_SRC_FILES})

if(MSVC)
	target_compile_options(${LIB_WFA} PUBLIC /W4 /WX /Zc:preprocessor /Wv:18)
endif()

if(CMAKE_CXX_COMPILER_ID STREQUAL "Clang")
	target_compile_options(${LIB_WFA} PUBLIC -Werror -W -pedantic -fexperimental-library )
	# set_target_properties(${LIB_SYNTH} PROPERTIES CXX_CLANG_TIDY "clang-tidy")
endif()

target_link_libraries(${LIB_WFA} PUBLIC
					fmt::fmt
					Kokkos::kokkos
)

set_property(TARGET ${LIB_WFA} PROPERTY CXX_STANDARD 20)

#
# WFA_TOOL
#

add_executable(wfa_tool "src/wfa_tool.cpp")
set_property(TARGET wfa_tool PROPERTY CXX_STANDARD 20)
target_link_libraries(wfa_tool PRIVATE ${LIB_WFA})

if(CMAKE_CXX_COMPILER_ID STREQUAL "Clang")
	target_compile_options(wfa_tool PUBLIC -Werror -W -pedantic -fexperimental-library )
	# set_target_properties(${LIB_SYNTH} PROPERTIES CXX_CLANG_TIDY "clang-tidy")
endif()

#
# TESTS
#

include(FetchContent)
# Fetch Catch2 library
FetchContent_Declare(
    Catch2
    GIT_REPOSITORY https://github.com/catchorg/Catch2.git
    GIT_TAG v3.7.1  # Replace with the latest version as needed
)
FetchContent_MakeAvailable(Catch2)

add_subdirectory(tests)

enable_testing()

add_executable(wfa_tests ${WFA_TEST_FILES})
target_link_libraries(wfa_tests PRIVATE Catch2::Catch2 ${LIB_WFA})

add_test(NAME NaiveAlignmentTests COMMAND naive_tests)