cmake_minimum_required(VERSION 3.20)
project(CCR)

set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
IF(CMAKE_SYSTEM_NAME MATCHES "Linux")
    set(CMAKE_CXX_COMPILER "/usr/bin/g++")
    set(CMAKE_C_COMPILER "/usr/bin/gcc")
ELSE()
    MESSAGE(FATAL_ERROR "You are not using a Linux system. Please change the compiler path manually in the CMakeLists.txt.")
ENDIF()
find_package(gflags)

add_definitions()

find_package(Threads REQUIRED)
SET(ROOT_DIR ${CMAKE_SOURCE_DIR})
SET(THIRD_DIR ${ROOT_DIR}/third)
SET(SPD_DIR ${THIRD_DIR}/spdlog/include)
SET(OREACH ${ROOT_DIR}/third/O-Reach)

# for oreach
include_directories(${OREACH})
include_directories(${OREACH}/app)
include_directories(${OREACH}/extern/argtable3-3.0.3)
include_directories(${OREACH}/extern/flat_hash_map)
include_directories(${OREACH}/lib)
include_directories(${OREACH}/extern/KaHIP/lib)
include_directories(${OREACH}/extern/KaHIP/lib/io)
include_directories(${OREACH}/extern/KaHIP/lib/tools)
include_directories(include ${SPD_DIR})

set(LIB_COMMON_SOURCE_FILES
  ${OREACH}/extern/KaHIP/lib/algorithms/strongly_connected_components.cpp
  ${OREACH}/extern/KaHIP/lib/io/graph_io.cpp
  ${OREACH}/extern/KaHIP/lib/io/graph_io.h
  ${OREACH}/extern/KaHIP/lib/tools/random_functions.cpp
  ${OREACH}/extern/KaHIP/lib/tools/macros_assertions.h
  ${OREACH}/extern/KaHIP/lib/tools/macros_common.h
  ${OREACH}/extern/KaHIP/lib/definitions.h
  ${OREACH}/extern/KaHIP/lib/data_structure/union_find.h
  ${OREACH}/app/query_generator.h
  ${OREACH}/lib/algorithms/topological_sort.cpp
  ${OREACH}/extern/argtable3-3.0.3/argtable3.c
  ${OREACH}/extern/argtable3-3.0.3/argtable3.h
  ${OREACH}/extern/flat_hash_map/flat_hash_map.hpp
  ${OREACH}/extern/flat_hash_map/bytell_hash_map.hpp
  ${OREACH}/extern/murmur_hash/MurmurHash3.cpp
  )
add_library(libcommon OBJECT ${LIB_COMMON_SOURCE_FILES})

set(LIBREACHABILITY_SOURCE_FILES OBJECT 
  ${OREACH}/lib/algorithms/bfs.cpp
  ${OREACH}/lib/reachability/oreach.cpp
  ${OREACH}/lib/reachability/matrix/full_reach.h
)
add_library(libreachability ${LIBREACHABILITY_SOURCE_FILES})

set(CMAKE_CXX_FLAGS "-O3 -fopenmp")
# set(CMAKE_CXX_FLAGS "-g -fopenmp")
set(CMAKE_CXX_STANDARD 14)

set(SOURCE_FILE Computation.cpp)
# aux_source_directory(${CMAKE_SOURCE_DIR}/src SRC)
file(GLOB_RECURSE SRC ${CMAKE_SOURCE_DIR}/src/*.cpp)


add_executable(Computation ${SOURCE_FILE} ${SRC} $<TARGET_OBJECTS:libcommon> $<TARGET_OBJECTS:libreachability>)
target_link_libraries (Computation gflags Threads::Threads)
add_subdirectory(tool)



file(GLOB_RECURSE ALL_SOURCE_FILES 
     ${CMAKE_SOURCE_DIR}/*.c 
     ${CMAKE_SOURCE_DIR}/*.cpp 
     ${CMAKE_SOURCE_DIR}/*.h)

# 排除 third 目录
list(FILTER ALL_SOURCE_FILES EXCLUDE REGEX "${CMAKE_SOURCE_DIR}/third/.*")

# 自定义命令以运行 clang-format
add_custom_command(OUTPUT clang_format
                   COMMAND clang-format
                   ARGS -i -style=file ${ALL_SOURCE_FILES}
                   COMMENT "Running clang-format on all source files except those under 'third/' directory"
                   VERBATIM)

# 创建伪目标
add_custom_target(format DEPENDS clang_format)