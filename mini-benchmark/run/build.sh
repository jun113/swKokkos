#!/usr/bin/bash

set -x
set -e

cd ..

CUR_DIR=$(pwd)

mkdir -p build

cd build

  /usr/sw/cmake-3.23.0-rc3-linux-x86_64/bin/cmake .. \
    -D CMAKE_C_COMPILER="swgcc" \
    -D CMAKE_CXX_COMPILER="swg++" \
    -D Kokkos_ENABLE_DEBUG_BOUNDS_CHECK=off \
    -D Kokkos_ENABLE_COMPILER_WARNINGS=off \
    -D Kokkos_ENABLE_DEBUG=Off \
    -D Kokkos_ENABLE_AGGRESSIVE_VECTORIZATION=on \
    -D Kokkos_ENABLE_SERIAL=On \
    -D Kokkos_ENABLE_ATHREAD=on \
    -D Kokkos_ATHREAD_FUNCTOR_DIR="/home/export/base/shisuan/swiap/online/wjl/swKokkos-paper/dev/Kokkos-Athread/src"


    #-D Kokkos_ENABLE_SERIAL=on
    # -D CMAKE_CXX_COMPILER="/home/wjl/workspace/test/20230501-tgrid_to_ugrid/kokkos-4.0.00/bin/nvcc_wrapper" \
    #-D CMAKE_CXX_COMPILER="${KOKKOS_SRC}/bin/nvcc_wrapper" \
    # -D Kokkos_ENABLE_CUDA_CONSTEXPR=on \
    # -D Kokkos_ENABLE_CUDA_LAMBDA=on \
    # -D Kokkos_ENABLE_CUDA=on \
    # -D Kokkos_ARCH_VOLTA70=on \

#make
#make VERBOSE=1 -j
make -j

cp swKokkosTest ../run/
