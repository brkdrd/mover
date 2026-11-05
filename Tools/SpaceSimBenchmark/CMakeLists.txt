cmake_minimum_required(VERSION 3.20 FATAL_ERROR)

# Record the exact compile commands of both translation units. This must be set
# before project()/the first target for the generators to honour it; the run
# script records the library and harness commands from it.
set(CMAKE_EXPORT_COMPILE_COMMANDS ON CACHE BOOL "Write compile_commands.json for provenance" FORCE)

# Standalone, portable subproject for the space-dynamics throughput benchmark.
#
# It builds the repository's Jolt library out of the existing Build/ project (the
# supported way to get the exact compiler/ISA definitions the engine is compiled
# with -- AVX2/SSE4.2/FMA/BMI etc. come from the Jolt target's PUBLIC interface)
# and links one small benchmark executable against it. Nothing outside this
# directory is modified.
#
# Usage (from a clean build directory):
#   cmake -S Tools/SpaceSimBenchmark -B <build-dir> -DCMAKE_BUILD_TYPE=Release
#   cmake --build <build-dir> -j2
#
# The options below pin the measurement configuration. They are set in the cache
# before Jolt's option() calls, so they win. A throughput run must not use the
# instrumenting profiler or the debug renderer (both change measured cost), and
# this is an optimized float-precision build with asserts off, as reported in the
# manifest.
#
# The generated SpaceSimBenchmarkBuildInfo.h records the compiler, the flags and
# definitions that were actually selected for the library and for the benchmark
# translation unit, the relevant cache values, and content fingerprints of the
# compiled harness inputs and of the Jolt sources linked into the binary. The
# library and the harness TU do NOT get identical flags (Build/CMakeLists.txt
# adds its flags in the Build/ directory scope, which does not reach this
# target), so both are recorded separately. compile_commands.json holds the exact
# commands; run_benchmark.sh records them next to the results.
project(SpaceSimBenchmark LANGUAGES CXX)

set(SPACESIMBENCH_JOLT_ROOT "${CMAKE_CURRENT_SOURCE_DIR}/../.." CACHE PATH "Repository root containing Jolt/ and Build/")

set(PROFILER_IN_DEBUG_AND_RELEASE OFF CACHE BOOL "Benchmark: no instrumenting profiler" FORCE)
set(PROFILER_IN_DISTRIBUTION OFF CACHE BOOL "Benchmark: no instrumenting profiler" FORCE)
set(DEBUG_RENDERER_IN_DEBUG_AND_RELEASE OFF CACHE BOOL "Benchmark: headless" FORCE)
set(DEBUG_RENDERER_IN_DISTRIBUTION OFF CACHE BOOL "Benchmark: headless" FORCE)
set(USE_ASSERTS OFF CACHE BOOL "Benchmark: Release without asserts" FORCE)
set(DOUBLE_PRECISION OFF CACHE BOOL "Benchmark: float position precision (engine default)" FORCE)
set(ENABLE_ALL_WARNINGS ON CACHE BOOL "Let the library build with its normal warning set" FORCE)

# Build only the physics library. As a subproject Build/CMakeLists.txt skips its
# sample/test/viewer targets automatically (it checks for being the root project).
add_subdirectory(${SPACESIMBENCH_JOLT_ROOT}/Build ${CMAKE_BINARY_DIR}/jolt)

add_executable(SpaceSimBenchmark src/main.cpp)
set_target_properties(SpaceSimBenchmark PROPERTIES
	CXX_STANDARD 17
	CXX_STANDARD_REQUIRED ON
	CXX_EXTENSIONS OFF)
if (MSVC)
	target_compile_options(SpaceSimBenchmark PRIVATE /W4)
else()
	target_compile_options(SpaceSimBenchmark PRIVATE -Wall -Wextra)
endif()
target_link_libraries(SpaceSimBenchmark PRIVATE Jolt)

# ---------------------------------------------------------------------------
# Build/provenance record embedded in the benchmark.
#
# Harness fingerprint: the workload source and the CMake file that selects the
# measurement configuration. Documentation, the run script and the tests are
# deliberately excluded so editing them does not invalidate a binary. The
# fingerprint is recomputed by run_benchmark.sh from the checked-out files and
# must match the binary's manifest.
# ---------------------------------------------------------------------------
set(_spacesim_fingerprint_sources src/main.cpp CMakeLists.txt)
set(_spacesim_source_manifest "")
foreach(_f IN LISTS _spacesim_fingerprint_sources)
	file(SHA256 "${CMAKE_CURRENT_SOURCE_DIR}/${_f}" _h)
	string(APPEND _spacesim_source_manifest "${_f}:${_h}\n")
endforeach()
string(SHA256 SPACESIMBENCH_SOURCE_SHA256 "${_spacesim_source_manifest}")

# Core fingerprint: content hash of every existing source file of the Jolt
# target. This is what the binary actually links, so it changes when the core
# changes even if the working-tree dirty count happens to look the same.
get_target_property(_jolt_sources Jolt SOURCES)
set(_spacesim_core_manifest "")
set(_spacesim_core_files 0)
foreach(_src IN LISTS _jolt_sources)
	if (NOT IS_ABSOLUTE "${_src}" OR NOT EXISTS "${_src}" OR IS_DIRECTORY "${_src}")
		continue()
	endif()
	file(RELATIVE_PATH _rel "${SPACESIMBENCH_JOLT_ROOT}" "${_src}")
	file(SHA256 "${_src}" _h)
	string(APPEND _spacesim_core_manifest "${_rel}:${_h}\n")
	math(EXPR _spacesim_core_files "${_spacesim_core_files} + 1")
endforeach()
string(SHA256 SPACESIMBENCH_CORE_SHA256 "${_spacesim_core_manifest}")

# Regenerate this header whenever a compiled input changes, so the embedded
# fingerprints cannot silently describe an older tree than the binary.
set(_spacesim_configure_depends "${CMAKE_CURRENT_SOURCE_DIR}/src/main.cpp" "${CMAKE_CURRENT_SOURCE_DIR}/CMakeLists.txt")
foreach(_src IN LISTS _jolt_sources)
	if (IS_ABSOLUTE "${_src}" AND EXISTS "${_src}" AND NOT IS_DIRECTORY "${_src}")
		list(APPEND _spacesim_configure_depends "${_src}")
	endif()
endforeach()
set_property(DIRECTORY APPEND PROPERTY CMAKE_CONFIGURE_DEPENDS ${_spacesim_configure_depends})

# Definitions propagated from the Jolt target (ISA selection, precision, ...).
# Entries that are generator expressions (e.g. NDEBUG for optimization configs)
# are not rendered here; the exact -D list is in the compile command recorded by
# run_benchmark.sh, and -DNDEBUG is visible in the flag strings below.
get_target_property(_jolt_defs Jolt INTERFACE_COMPILE_DEFINITIONS)
if (NOT _jolt_defs)
	set(_jolt_defs "")
endif()
set(_spacesim_defs_text "")
foreach(_d IN LISTS _jolt_defs)
	if (NOT _d MATCHES "^[A-Za-z_][A-Za-z0-9_]*$")
		continue()
	endif()
	string(APPEND _spacesim_defs_text "-D${_d} ")
endforeach()
string(STRIP "${_spacesim_defs_text}" _spacesim_defs_text)

# Collapse repeated flags (CMake lists the PUBLIC Jolt options in both
# COMPILE_OPTIONS and INTERFACE_COMPILE_OPTIONS) while preserving order.
function(_spacesim_dedupe_words in_text out_var)
	string(REPLACE " " ";" _tokens "${in_text}")
	set(_unique "")
	foreach(_t IN LISTS _tokens)
		if (_t STREQUAL "")
			continue()
		endif()
		list(FIND _unique "${_t}" _index)
		if (_index EQUAL -1)
			list(APPEND _unique "${_t}")
		endif()
	endforeach()
	string(REPLACE ";" " " _deduped "${_unique}")
	set(${out_var} "${_deduped}" PARENT_SCOPE)
endfunction()

# Flags as configured for each translation unit. Build/CMakeLists.txt sets its
# flags on the CMAKE_CXX_FLAGS variable of the Build/ directory scope, so the
# library gets them and this target does not; both sets are recorded.
get_directory_property(_jolt_dir_cxx_flags DIRECTORY "${SPACESIMBENCH_JOLT_ROOT}/Build" DEFINITION CMAKE_CXX_FLAGS)
get_target_property(_jolt_compile_options Jolt COMPILE_OPTIONS)
get_target_property(_jolt_interface_options Jolt INTERFACE_COMPILE_OPTIONS)
get_target_property(_bench_compile_options SpaceSimBenchmark COMPILE_OPTIONS)
if (NOT _jolt_dir_cxx_flags)
	set(_jolt_dir_cxx_flags "")
endif()
if (NOT _jolt_compile_options)
	set(_jolt_compile_options "")
endif()
if (NOT _jolt_interface_options)
	set(_jolt_interface_options "")
endif()
if (NOT _bench_compile_options)
	set(_bench_compile_options "")
endif()

if (CMAKE_BUILD_TYPE)
	string(TOUPPER "${CMAKE_BUILD_TYPE}" _cfg)
else()
	set(_cfg "RELEASE")
endif()
set(_config_cxx_flags "${CMAKE_CXX_FLAGS_${_cfg}}")

string(REPLACE ";" " " _jolt_dir_cxx_flags "${_jolt_dir_cxx_flags}")
string(REPLACE ";" " " _jolt_compile_options "${_jolt_compile_options}")
string(REPLACE ";" " " _jolt_interface_options "${_jolt_interface_options}")
string(REPLACE ";" " " _bench_compile_options "${_bench_compile_options}")

set(_bench_flags "${CMAKE_CXX_FLAGS} ${_config_cxx_flags} -std=c++17 ${_jolt_interface_options} ${_bench_compile_options}")
set(_lib_flags "${CMAKE_CXX_FLAGS} ${_config_cxx_flags} ${_jolt_dir_cxx_flags} -std=c++17 ${_jolt_compile_options}")
string(STRIP "${_bench_flags}" _bench_flags)
string(STRIP "${_lib_flags}" _lib_flags)
string(REGEX REPLACE "[ ]+" " " _bench_flags "${_bench_flags}")
string(REGEX REPLACE "[ ]+" " " _lib_flags "${_lib_flags}")
_spacesim_dedupe_words("${_bench_flags}" _bench_flags)
_spacesim_dedupe_words("${_lib_flags}" _lib_flags)

# Facts that are easy to get wrong when the parent presets are assumed.
if (_bench_flags MATCHES "-ffast-math|/fp:fast" OR _lib_flags MATCHES "-ffast-math|/fp:fast")
	set(_fast_math "true")
else()
	set(_fast_math "false")
endif()
if (INTERPROCEDURAL_OPTIMIZATION)
	set(_lib_lto "on")
else()
	set(_lib_lto "off")
endif()
set(_bench_lto "off") # the harness TU is not compiled with LTO; see _bench_flags
if (CROSS_PLATFORM_DETERMINISTIC)
	set(_cross_platform_deterministic "on")
else()
	set(_cross_platform_deterministic "off")
endif()

set(_cache_record "")
foreach(_c INTERPROCEDURAL_OPTIMIZATION CROSS_PLATFORM_DETERMINISTIC DOUBLE_PRECISION USE_ASSERTS
		USE_AVX2 USE_FMADD USE_SSE4_2 PROFILER_IN_DEBUG_AND_RELEASE DEBUG_RENDERER_IN_DEBUG_AND_RELEASE
		ENABLE_ALL_WARNINGS JPH_BUILD_SHARED_LIBS)
	if (DEFINED ${_c})
		string(APPEND _cache_record "${_c}=${${_c}};")
	endif()
endforeach()
string(STRIP "${_cache_record}" _cache_record)

# Configure-time git state. The dirty-entry count is not an identity of the tree
# (it does not say which lines changed), which is why the content fingerprints
# above exist; both are recorded so the scope is visible.
set(_git_head "unknown")
set(_git_dirty_entries "unknown")
find_package(Git QUIET)
if (GIT_FOUND)
	execute_process(COMMAND ${GIT_EXECUTABLE} -C "${SPACESIMBENCH_JOLT_ROOT}" rev-parse HEAD
		OUTPUT_VARIABLE _git_head OUTPUT_STRIP_TRAILING_WHITESPACE ERROR_QUIET RESULT_VARIABLE _git_rc)
	if (NOT _git_rc EQUAL 0)
		set(_git_head "unknown")
	endif()
	execute_process(COMMAND ${GIT_EXECUTABLE} -C "${SPACESIMBENCH_JOLT_ROOT}" status --porcelain
		OUTPUT_VARIABLE _git_status OUTPUT_STRIP_TRAILING_WHITESPACE ERROR_QUIET RESULT_VARIABLE _git_rc)
	if (_git_rc EQUAL 0)
		string(REGEX MATCHALL "[^\n]+" _git_status_lines "${_git_status}")
		list(LENGTH _git_status_lines _git_dirty_entries)
	else()
		set(_git_dirty_entries "unknown")
	endif()
endif()

set(_spacesim_build_manifest "${_spacesim_source_manifest}")
string(APPEND _spacesim_build_manifest "core:${SPACESIMBENCH_CORE_SHA256}\n")
string(APPEND _spacesim_build_manifest "defs:${_spacesim_defs_text}\nbench:${_bench_flags}\nlib:${_lib_flags}\ncache:${_cache_record}\n")
string(SHA256 SPACESIMBENCH_BUILD_SHA256 "${_spacesim_build_manifest}")

file(WRITE "${CMAKE_CURRENT_BINARY_DIR}/SpaceSimBenchmarkBuildInfo.h"
"// Generated by Tools/SpaceSimBenchmark/CMakeLists.txt -- do not edit.
#pragma once
#define SPACESIMBENCH_COMPILER \"${CMAKE_CXX_COMPILER_ID} ${CMAKE_CXX_COMPILER_VERSION}\"
#define SPACESIMBENCH_BUILD_TYPE \"${CMAKE_BUILD_TYPE}\"
#define SPACESIMBENCH_SOURCE_SHA256 \"${SPACESIMBENCH_SOURCE_SHA256}\"
#define SPACESIMBENCH_CORE_SHA256 \"${SPACESIMBENCH_CORE_SHA256}\"
#define SPACESIMBENCH_BUILD_SHA256 \"${SPACESIMBENCH_BUILD_SHA256}\"
#define SPACESIMBENCH_CORE_SOURCE_FILES ${_spacesim_core_files}
#define SPACESIMBENCH_GIT_HEAD \"${_git_head}\"
#define SPACESIMBENCH_GIT_DIRTY_ENTRIES ${_git_dirty_entries}
#define SPACESIMBENCH_DEFINITIONS \"${_spacesim_defs_text}\"
#define SPACESIMBENCH_BENCH_COMPILE_FLAGS \"${_bench_flags}\"
#define SPACESIMBENCH_LIB_COMPILE_FLAGS \"${_lib_flags}\"
#define SPACESIMBENCH_CACHE \"${_cache_record}\"
#define SPACESIMBENCH_LIB_LTO \"${_lib_lto}\"
#define SPACESIMBENCH_BENCH_LTO \"${_bench_lto}\"
#define SPACESIMBENCH_CROSS_PLATFORM_DETERMINISTIC \"${_cross_platform_deterministic}\"
#define SPACESIMBENCH_FAST_MATH \"${_fast_math}\"
")
target_include_directories(SpaceSimBenchmark PRIVATE "${CMAKE_CURRENT_BINARY_DIR}")

# ---------------------------------------------------------------------------
# Regression tests: CLI validation and artifact I/O failure handling, driven by
# the real binary as a subprocess. Run with `ctest` or
#   python3 Tools/SpaceSimBenchmark/tests/test_cli_and_artifacts.py <binary>
# ---------------------------------------------------------------------------
enable_testing()
find_package(Python3 COMPONENTS Interpreter QUIET)
if (Python3_Interpreter_FOUND)
	add_test(NAME SpaceSimBenchmark.CliAndArtifacts
		COMMAND ${Python3_EXECUTABLE} "${CMAKE_CURRENT_SOURCE_DIR}/tests/test_cli_and_artifacts.py" "$<TARGET_FILE:SpaceSimBenchmark>")
	set_tests_properties(SpaceSimBenchmark.CliAndArtifacts PROPERTIES TIMEOUT 600)
endif()
