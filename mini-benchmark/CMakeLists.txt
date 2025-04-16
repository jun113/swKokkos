#CMake >= 3.18: Fortran linkage. This does not affect most mixed Fortran/Kokkos builds.
#New in version 3.8: Added CUDA support.
#In CMake 3.12 and above the find_package() 
#command now searches prefixes specified 
#by the _ROOT CMake variable and the ROOT environment variable. 
# Taishan
#cmake_minimum_required (VERSION 3.16)
cmake_minimum_required (VERSION 3.18)

#project name
set(PROJECT_NAME "swKokkosTest")
project(${PROJECT_NAME})

enable_language(C CXX)

# Debug or Release
#set(CMAKE_BUILD_TYPE               Debug)
set(CMAKE_BUILD_TYPE               Release)
# set(EXECUTABLE_OUTPUT_PATH         ${CMAKE_SOURCE_DIR}/build/bin)

#set(CMAKE_CXX_COMPILER "nvcc")
#set(CMAKE_CXX_FLAGS	"-O3 -Wall -gencode=arch=compute_70,code=sm_70")
# set(CMAKE_CXX_FLAGS	"-O3 -Wall")
set(CMAKE_C_FLAGS	"-O3")
set(CMAKE_CXX_FLAGS	"-O3 -fwrapv")

set(CMAKE_CXX_STANDARD 17)

set(CMAKE_EXE_LINKER_FLAGS	"${CMAKE_EXE_LINKER_FLAGS} -rdynamic -finstrument-functions")

# GPTL
#include_directories(/home/wjl/software/GPTL/8.1.1/include)
#link_directories(/home/wjl/software/GPTL/8.1.1/lib)

# CUDA
#include_directories(/usr/local/cuda/include)
#link_directories(/usr/local/cuda/lib64)

# Athread
include_directories(/usr/sw/swgcc/swgcc710-tools-SEA/shared_include)
link_directories(/usr/sw/swgcc/swgcc710-tools-SEA-1307/usr/sw_64sw6a-sunway-linux-gnu/lib)

file(GLOB SRC-HOST  LIST_DIRECTORIES false RELATIVE ${CMAKE_SOURCE_DIR} 
		src/main.cpp 
		src/athread-impl/athread_axpy.c
		src/athread-impl/athread_3d_mult.c
		src/athread-impl/athread_dot_product_2d.c
		src/athread-impl/athread_stencil.c
		src/athread-impl/athread_matrix_add.c)
file(GLOB SRC-SLAVE LIST_DIRECTORIES false RELATIVE ${CMAKE_SOURCE_DIR} 
		src/athread-impl/athread_axpy_slave.c
		src/athread-impl/athread_3d_mult_slave.c
		src/athread-impl/athread_dot_product_2d_slave.c
		src/athread-impl/athread_stencil_slave.c
		src/athread-impl/athread_matrix_add_slave.c
		)

set_source_files_properties(${SRC-HOST} PROPERTIES COMPILE_FLAGS "-mhost -fPIC")
set_source_files_properties(${SRC-SLAVE} PROPERTIES COMPILE_FLAGS "-mslave -fPIC -DKOKKOS_ENABLE_ATHREAD_SLAVE")

set(CMAKE_EXE_LINKER_FLAGS	"${CMAKE_EXE_LINKER_FLAGS} -mhybrid -fPIC")

add_subdirectory(${CMAKE_SOURCE_DIR}/Kokkos/swKokkos-4.1.00)

add_executable(${PROJECT_NAME} ${SRC-HOST} ${SRC-SLAVE})

#target_link_libraries(${PROJECT_NAME} gptl libmpifort.so.12 nvToolsExt  Kokkos::kokkos)
target_link_libraries(${PROJECT_NAME} Kokkos::kokkos  m_slave m)
#  mfma
