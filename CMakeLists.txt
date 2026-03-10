# Set minimum required CMake version
cmake_minimum_required(VERSION 3.16)

# Append the path to the custom cmake modules
list(APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake/modules)

# Start project, set its name to shared-memory-hypergraph-mincut and set languages
project(shared-memory-hypergraph-mincut C CXX)

# Force compiler to run in C++14 mode
set(CMAKE_CXX_STANDARD 14)
# Require the CXX_STANDARD
set(CMAKE_CXX_STANDARD_REQUIRED ON)
# Disable compiler specific extensions
set(CMAKE_CXX_EXTENSIONS OFF)

# if no build mode is specified build in release mode
if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  set(CMAKE_BUILD_TYPE "Release")
endif()

# Find TBB library (taken from Mt-KaHyPar)
# NB: This is only necessary if you want to import the Mt-KaHyPar headers directly
find_package(TBB REQUIRED)
if (TBB_FOUND)
  if ((TBB_INTERFACE_VERSION GREATER_EQUAL 12050) OR (NOT KAHYPAR_ENFORCE_MINIMUM_TBB_VERSION) OR KAHYPAR_CI_BUILD)
    include_directories(SYSTEM ${TBB_INCLUDE_DIRS} ${TBB_MALLOC_INCLUDE_DIRS})
    link_libraries(${TBB_LIBRARIES} ${TBB_MALLOC_LIBRARIES})
    MESSAGE(STATUS "TBB Version: ${TBB_VERSION_MAJOR}.${TBB_VERSION_MINOR} (Interface Version=${TBB_INTERFACE_VERSION})")
    MESSAGE(STATUS "Found TBB library: inc=${TBB_INCLUDE_DIRS}, lib=${TBB_LIBRARIES}")
    MESSAGE(STATUS "Found TBB Malloc Proxy library: inc=${TBB_MALLOC_PROXY_INCLUDE_DIRS}, lib=${TBB_MALLOC_PROXY_LIBRARIES}")
  ELSE ()
    MESSAGE(FATAL_ERROR "
      Your current TBB version is too old. Your interface version is ${TBB_INTERFACE_VERSION}, 
      but the minimum required version is 12050 (OneTBB 2021.5).")
  ENDIF()
ELSE ()
  MESSAGE(FATAL_ERROR "TBB library not found. Install TBB on your system.")
ENDIF()

# Find Gurobi library (taken from KaHIP)
find_package(Gurobi REQUIRED)
if(GUROBI_FOUND)
  MESSAGE("Using Gurobi for ILP solver in ilp_improve")
  include_directories(SYSTEM ${GUROBI_INCLUDE_DIR})
  link_libraries(${GUROBI_LIBRARIES})
ELSE()
MESSAGE(FATAL_ERROR "Gurobui library not found. Install Gurobui on your system.")
ENDIF()

# Find Boost Program Options Library (taken from Mt-KaHyPar)
find_package(Boost 1.69 REQUIRED COMPONENTS program_options)
if(Boost_FOUND)
  include_directories(SYSTEM ${Boost_INCLUDE_DIRS})
  set(KAHYPAR_INCLUDE_DIRS ${Boost_INCLUDE_DIRS} ${KAHYPAR_INCLUDE_DIRS})
  set(KAHYPAR_LINK_LIBRARIES ${Boost_LIBRARIES} ${KAHYPAR_LINK_LIBRARIES})
  message(STATUS "Boost Include: ${Boost_INCLUDE_DIRS}")
  message(STATUS "Boost Library Dirs: ${Boost_LIBRARY_DIRS}")
  message(STATUS "Boost Libraries: ${Boost_LIBRARIES}")
  if(WIN32)
    add_definitions(-DBOOST_ALL_NO_LIB)
    add_definitions(-DBOOST_PROGRAM_OPTIONS_DYN_LINK=1)
  endif()
else()
  MESSAGE(FATAL_ERROR "
    Boost not found. Install Boost on your system.")
endif()

# Find the Threads library
find_package(Threads REQUIRED)
if(Threads_FOUND)
  message(STATUS "Found Threads: ${CMAKE_THREAD_LIBS_INIT}")
else()
  message(FATAL_ERROR "Threads library not found. Install Threads on your system.")
endif()

# Include directories
include_directories(${CMAKE_CURRENT_SOURCE_DIR})
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/extern)
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/extern/kahip)
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/extern/mt-kahypar)
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/extern/mt-kahypar-library)
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/extern/kahypar-shared-resources)
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/extern/argtable3-3.2.2)

# Create the dependencies
set(DEP_common  ${CMAKE_CURRENT_SOURCE_DIR}/lib/parse_parameters/parse_parameters.cpp
                ${CMAKE_CURRENT_SOURCE_DIR}/extern/argtable3-3.2.2/argtable3.c
                ${CMAKE_CURRENT_SOURCE_DIR}/lib/utils/random.cpp)

set(DEP_kernelizer  ${DEP_common}
                    ${CMAKE_CURRENT_SOURCE_DIR}/app/mincut_kernelizer.cpp
                    ${CMAKE_CURRENT_SOURCE_DIR}/lib/solvers/kernelizer.cpp
                    ${CMAKE_CURRENT_SOURCE_DIR}/lib/coarsening/label_propagation.cpp
                    ${CMAKE_CURRENT_SOURCE_DIR}/lib/solvers/ilp.cpp
                    ${CMAKE_CURRENT_SOURCE_DIR}/lib/solvers/submodular.cpp)

set(DEP_ilp ${DEP_common}
            ${CMAKE_CURRENT_SOURCE_DIR}/app/mincut_ilp.cpp
            ${CMAKE_CURRENT_SOURCE_DIR}/lib/solvers/ilp.cpp)

set(DEP_trimmer  ${DEP_common}
                    ${CMAKE_CURRENT_SOURCE_DIR}/app/mincut_trimmer.cpp
                    ${CMAKE_CURRENT_SOURCE_DIR}/lib/solvers/submodular.cpp
                    ${CMAKE_CURRENT_SOURCE_DIR}/lib/trimmer/trimmer.cpp)

set(DEP_submodular  ${DEP_common}
                    ${CMAKE_CURRENT_SOURCE_DIR}/app/mincut_submodular.cpp
                    ${CMAKE_CURRENT_SOURCE_DIR}/lib/solvers/submodular.cpp)

set(DEP_dumbbell_generator  ${DEP_common}
                            ${CMAKE_CURRENT_SOURCE_DIR}/app/generators/dumbbell_generator.cpp)

set(DEP_kcore_generator ${DEP_common}
                        ${CMAKE_CURRENT_SOURCE_DIR}/lib/decomposition/core_decomposition.cpp
                        ${CMAKE_CURRENT_SOURCE_DIR}/app/generators/kcore_generator.cpp
                        ${CMAKE_CURRENT_SOURCE_DIR}/lib/solvers/kernelizer.cpp
                        ${CMAKE_CURRENT_SOURCE_DIR}/lib/coarsening/label_propagation.cpp
                        ${CMAKE_CURRENT_SOURCE_DIR}/lib/solvers/ilp.cpp
                        ${CMAKE_CURRENT_SOURCE_DIR}/lib/solvers/submodular.cpp)

# Create executable programs
add_executable(heicut_kernelizer           ${DEP_kernelizer})
add_executable(heicut_kernelizer_parallel  ${DEP_kernelizer})
add_executable(heicut_ilp                  ${DEP_ilp})
add_executable(heicut_ilp_parallel         ${DEP_ilp})
add_executable(heicut_trimmer              ${DEP_trimmer})
add_executable(heicut_submodular           ${DEP_submodular})
add_executable(heicut_submodular_parallel  ${DEP_submodular})
add_executable(heicut_dumbbell_generator   ${DEP_dumbbell_generator})
add_executable(heicut_kcore_generator      ${DEP_kcore_generator})

# Link the Mt-KaHyPar library
target_link_libraries(heicut_kernelizer          ${CMAKE_CURRENT_SOURCE_DIR}/extern/mt-kahypar-library/libmtkahypar.so ${Boost_LIBRARIES})
target_link_libraries(heicut_kernelizer_parallel ${CMAKE_CURRENT_SOURCE_DIR}/extern/mt-kahypar-library/libmtkahypar.so ${Boost_LIBRARIES} ${CMAKE_THREAD_LIBS_INIT})
target_link_libraries(heicut_ilp                 ${CMAKE_CURRENT_SOURCE_DIR}/extern/mt-kahypar-library/libmtkahypar.so ${Boost_LIBRARIES})
target_link_libraries(heicut_ilp_parallel        ${CMAKE_CURRENT_SOURCE_DIR}/extern/mt-kahypar-library/libmtkahypar.so ${Boost_LIBRARIES})
target_link_libraries(heicut_trimmer             ${CMAKE_CURRENT_SOURCE_DIR}/extern/mt-kahypar-library/libmtkahypar.so ${Boost_LIBRARIES})
target_link_libraries(heicut_submodular          ${CMAKE_CURRENT_SOURCE_DIR}/extern/mt-kahypar-library/libmtkahypar.so ${Boost_LIBRARIES})
target_link_libraries(heicut_submodular_parallel ${CMAKE_CURRENT_SOURCE_DIR}/extern/mt-kahypar-library/libmtkahypar.so ${Boost_LIBRARIES} ${CMAKE_THREAD_LIBS_INIT})
target_link_libraries(heicut_kcore_generator     ${CMAKE_CURRENT_SOURCE_DIR}/extern/mt-kahypar-library/libmtkahypar.so ${Boost_LIBRARIES})

# Add the parallelization flag
target_compile_definitions(heicut_ilp_parallel        PUBLIC SMHM_PARALLEL)
target_compile_definitions(heicut_kernelizer_parallel PUBLIC SMHM_PARALLEL)
target_compile_definitions(heicut_submodular_parallel PUBLIC SMHM_PARALLEL)