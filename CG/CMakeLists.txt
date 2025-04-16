cmake_minimum_required (VERSION 3.18)

set(PROJECT_NAME			"CG")
project (${PROJECT_NAME})
#find_package(Kokkos REQUIRED)

enable_language(C CXX)

# Debug or Release
set(CMAKE_BUILD_TYPE			          Debug)
set(EXECUTABLE_OUTPUT_PATH		      ${CMAKE_SOURCE_DIR}/build/bin)
set(LIBRARY_OUTPUT_PATH			        ${CMAKE_SOURCE_DIR}/build/bin)
#set(CMAKE_Fortran_MODULE_DIRECTORY	${CMAKE_SOURCE_DIR}/build/mod)


set(CMAKE_CXX_STANDARD 17)

if(${MACHINE} STREQUAL "419")
    set(CMAKE_CXX_FLAGS         "-O3 --gpu-max-threads-per-block=1024 -Wunused-result")
elseif(${MACHINE} STREQUAL "OCEANLIGHT")
    set(CMAKE_C_FLAGS           "-O3 -mftz -mieee -faddress_align=64")
    set(CMAKE_CXX_FLAGS         "-O3 -mftz -mieee -faddress_align=64 -std=gnu++1z")
else ()
    set(CMAKE_CXX_FLAGS         "-O3 -g")
endif ()

include_directories(${CMAKE_SOURCE_DIR}/src/inc)

if(${MACHINE} STREQUAL "IAP_GPU")

# netcdf
    link_directories(/usr/local/netcdf/lib)
    include_directories(/usr/local/netcdf/include)
# MPI
    link_directories(/home/wjl/software/oneAPI/mpi/2021.2.0/lib)
    include_directories(/home/wjl/software/oneAPI/mpi/2021.2.0/include)

# PAPI
    include_directories(/home/wjl/software/PAPI/7.0.1/include)
    link_directories(/home/wjl/software/PAPI/7.0.1/lib)
# GPTL
    include_directories(/home/wjl/software/GPTL/8.1.1/include)
    link_directories(/home/wjl/software/GPTL/8.1.1/lib)
# CUDA
    include_directories(/usr/local/cuda/include)
    link_directories(/usr/local/cuda/lib64)

elseif(${MACHINE} STREQUAL "419")
# MPI
    include_directories(/opt/hpc/software/mpi/hpcx/v2.4.1/intel-2017.5.239/include)
    link_directories(/opt/hpc/software/mpi/hpcx/v2.4.1/intel-2017.5.239/lib)
# GPTL
    include_directories(/public/home/esm7/software/GPTL/8.1.1/include)
    link_directories(/public/home/esm7/software/GPTL/8.1.1/lib)

elseif(${MACHINE} STREQUAL "KUNSHAN")

elseif(${MACHINE} STREQUAL "TAISHAN")
# Fortran Compiler
    link_directories(/workspace/public/software/compilers/gcc/kunpenggcc/9.3.1/gcc-9.3.1-2021.03-aarch64-linux/lib)
    link_directories(/workspace/public/software/compilers/gcc/kunpenggcc/9.3.1/gcc-9.3.1-2021.03-aarch64-linux/lib64)

# MPI
    include_directories(/workspace/public/software/mpi/openmpi/4.0.1_kunpenggcc9.3.1/include)
    link_directories(/workspace/public/software/mpi/openmpi/4.0.1_kunpenggcc9.3.1/lib)

# NETCDF
    include_directories(/workspace/home/migration/chengqian/huangjian/Build_WRF/library/netcdf/include)
    link_directories(/workspace/home/migration/chengqian/huangjian/Build_WRF/library/netcdf/lib)

    include_directories(/workspace/public/software/libs/netcdf/netcdf_c_4.7.4_fortran_4.5.3/kunpenggcc9.3.1_hmpi1.1.1/include)
    link_directories(/workspace/public/software/libs/netcdf/netcdf_c_4.7.4_fortran_4.5.3/kunpenggcc9.3.1_hmpi1.1.1/lib)
    
# HDF5
    include_directories(/workspace/home/migration/chengqian/huangjian/Build_WRF/library/hdf5/include)
    link_directories(/workspace/home/migration/chengqian/huangjian/Build_WRF/library/hdf5/lib)

    include_directories(/workspace/public/software/libs/hdf5/1.12.0/kunpenggcc9.3.1_hmpi1.1.1/include)
    link_directories(/workspace/public/software/libs/hdf5/1.12.0/kunpenggcc9.3.1_hmpi1.1.1/lib)
# GPTL
    include_directories(/workspace/home/migration/chengqian/wjl/GPTL/8.1.1/include)
    link_directories(/workspace/home/migration/chengqian/wjl/GPTL/8.1.1/lib)
# PAPI
    include_directories(/workspace/home/migration/chengqian/wjl/PAPI/7.0.1/include)
    link_directories(/workspace/home/migration/chengqian/wjl/PAPI/7.0.1/lib)
elseif(${MACHINE} STREQUAL "OCEANLIGHT")

    # NetCDF
    include_directories(/usr/sw/yyzlib/local-saca/include)
    link_directories(/usr/sw/yyzlib/local-saca/lib)

    # MPI
    include_directories(/usr/sw/mpi/mpi_current/include)

    # Athread
    include_directories(/usr/sw/swgcc/swgcc710-tools-SEA/shared_include)
    link_directories(/usr/sw/swgcc/swgcc710-tools-SEA-1307/usr/sw_64sw6a-sunway-linux-gnu/lib)

    # include_directories(/usr/sw/penv/include)
    # link_directories(/usr/sw/penv/lib)

    # include_directories(/usr/sw/yyzlib/local-saca/include)
    # link_directories(/usr/sw/yyzlib/local-saca/lib)
    
endif()

file(GLOB SRC-PCG	LIST_DIRECTORIES false RELATIVE ${CMAKE_SOURCE_DIR} 
    src/main.cpp
    src/pcg.cpp
    )

#file(GLOB SRC-PCG-HOST LIST_DIRECTORIES false RELATIVE ${CMAKE_SOURCE_DIR}
#    source/athread-impl/athread_lbm_host.c)

file(GLOB SRC-PCG-SLAVE LIST_DIRECTORIES false RELATIVE ${CMAKE_SOURCE_DIR}
    src/athread-impl/slave.c)

set_source_files_properties(${SRC-PCG} 
        PROPERTIES COMPILE_FLAGS "-mhost -fPIC -fpermissive")
set_source_files_properties(${SRC-PCG-HOST} 
        PROPERTIES COMPILE_FLAGS "-mhost -fPIC -fpermissive")
set_source_files_properties(${SRC-PCG-SLAVE} 
        PROPERTIES COMPILE_FLAGS "-mslave -fPIC -msimd")

set(CMAKE_EXE_LINKER_FLAGS	"${CMAKE_EXE_LINKER_FLAGS} -rdynamic -finstrument-functions")
set(CMAKE_EXE_LINKER_FLAGS	"${CMAKE_EXE_LINKER_FLAGS} -mhybrid -fPIC")

add_subdirectory(${CMAKE_SOURCE_DIR}/swKokkos-4.1.00)
#add_subdirectory(${CMAKE_SOURCE_DIR}/kokkos-4.5.01)

add_executable(${PROJECT_NAME} 
    ${SRC-PCG} 
    ${SRC-PCG-HOST} 
    ${SRC-PCG-SLAVE})

#set Fortran compiler as the linker
# set_target_properties(${PROJECT_NAME} PROPERTIES LINKER_LANGUAGE Fortran)

if(${MACHINE} STREQUAL "IAP_GPU")
    # target_link_libraries(${PROJECT_NAME} nvToolsExt)
elseif(${MACHINE} STREQUAL "OCEANLIGHT")
    target_link_libraries(${PROJECT_NAME} hdf5_hl hdf5 dl z m_slave m)
    set(CMAKE_C_IMPLICIT_LINK_LIBRARIES "stdc++")
    set(CMAKE_CXX_IMPLICIT_LINK_LIBRARIES "stdc++")
    set(CMAKE_C_IMPLICIT_LINK_LIBRARIES "-Wl,--as-needed -Wl,-Bdynamic")
endif()
    #target_link_libraries(${PROJECT_NAME} hdf5_hl hdf5 dl z m_slave m -lrt)

target_link_libraries(${PROJECT_NAME} Kokkos::kokkos)

# configure_file(${CMAKE_SOURCE_DIR}/run/def-undef.h.in     ${CMAKE_SOURCE_DIR}/src/head/def-undef.h     @ONLY)
# configure_file(${CMAKE_SOURCE_DIR}/run/kokkos_config.h.in ${CMAKE_SOURCE_DIR}/src/head/kokkos_config.h @ONLY)

if(${MACHINE} STREQUAL "IAP_GPU")
configure_file(${CMAKE_SOURCE_DIR}/script/run-iap-gpu.in    ${CMAKE_SOURCE_DIR}/script/run       @ONLY)
elseif(${MACHINE} STREQUAL "419")
configure_file(${CMAKE_SOURCE_DIR}/script/run-419.in        ${CMAKE_SOURCE_DIR}/script/run       @ONLY)
configure_file(${CMAKE_SOURCE_DIR}/script/dcu-bind-419.in   ${CMAKE_SOURCE_DIR}/script/dcu_bind  @ONLY)
configure_file(${CMAKE_SOURCE_DIR}/script/submit-419.in     ${CMAKE_SOURCE_DIR}/script/submit    @ONLY)
elseif(${MACHINE} STREQUAL "KUNSHAN")
#TODO
elseif(${MACHINE} STREQUAL "TAISHAN")
configure_file(${CMAKE_SOURCE_DIR}/script/run-taishan.in    ${CMAKE_SOURCE_DIR}/script/run        @ONLY)
configure_file(${CMAKE_SOURCE_DIR}/script/submit-taishan.in ${CMAKE_SOURCE_DIR}/script/submit     @ONLY)
elseif(${MACHINE} STREQUAL "OCEANLIGHT")
configure_file(${CMAKE_SOURCE_DIR}/script/run-ocean-light.in ${CMAKE_SOURCE_DIR}/script/run       @ONLY)
configure_file(${CMAKE_SOURCE_DIR}/script/submit-ocean-light.in ${CMAKE_SOURCE_DIR}/script/submit @ONLY)
endif()
