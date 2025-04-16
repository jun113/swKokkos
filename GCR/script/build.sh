#!/usr/bin/bash

#419 KUNSHAN IAP_GPU TAISHAN OCEANLIGHT
MACHINE="IAP_GPU"
MACHINE="OCEANLIGHT"

# CUR_DATE=$(date +'%Y%m%d-%H%M%S')

cd ../build

if [ "${MACHINE}" == "IAP_GPU" ]; then
	echo "IAP GPU"
	cmake	.. \
		-D Kokkos_ENABLE_SERIAL=On \
		-D Kokkos_ENABLE_CUDA=off \
    -D Kokkos_ARCH_VOLTA70=off \
    -D Kokkos_ENABLE_CUDA_CONSTEXPR=off \
		-D Kokkos_ENABLE_COMPILER_WARNINGS=off  \
		-D CMAKE_CXX_COMPILER="g++" \
		-D CMAKE_C_COMPILER="gcc" \
		-D CMAKE_Fortran_COMPILER="mpiifort" \
		-D MACHINE=${MACHINE}

		-D CMAKE_CXX_COMPILER="/home/wjl/workspace/swKokkos/lbm3d/lbm3d-kokkos/lbm3d/kokkos-4.5.01/bin/nvcc_wrapper" \
    # -D Kokkos_ENABLE_CUDA_LAMBDA=off \
		# -D CMAKE_CXX_COMPILER="/home/wjl/workspace/module/density/20231128/kokkos-4.1.00-athread/bin/nvcc_wrapper" \
elif [ "${MACHINE}" == "419" ]; then
	echo "419"
	module purge
	module load compiler/intel/2017.5.239
	module load mpi/hpcx/2.7.4/intel-2017.5.239
	module load compiler/rocm/dtk/22.10.1
	module list

	cmake	.. \
		-D Kokkos_ENABLE_SERIAL=On \
		-D Kokkos_ENABLE_HIP=off \
		-D Kokkos_ARCH_VEGA906=Off \
		-D Kokkos_ENABLE_COMPILER_WARNINGS=off  \
		-D CMAKE_C_COMPILER="hipcc" \
		-D CMAKE_CXX_COMPILER="hipcc" \
		-D CMAKE_Fortran_COMPILER="mpif90" \
		-D LOG_DIR=${LOG_DIR} \
		-D LICOM_ENABLE_FORTRAN_COMPILER=${LICOM_ENABLE_FORTRAN_COMPILER} \
    -D LICOM_ENABLE_VERSION=${LICOM_ENABLE_VERSION} \
		-D MACHINE=${MACHINE}

elif [ "${MACHINE}" == "TAISHAN" ]; then
	echo "Tai Shan"

elif [ "${MACHINE}" == "OCEANLIGHT" ]; then
	echo "Ocean Light"
  /usr/sw/cmake-3.23.0-rc3-linux-x86_64/bin/cmake .. \
    -D CMAKE_C_COMPILER="swgcc" \
    -D CMAKE_CXX_COMPILER="mpic++" \
    -D CMAKE_Fortran_COMPILER="mpif90" \
    -D Kokkos_ENABLE_DEBUG_BOUNDS_CHECK=off \
    -D Kokkos_ENABLE_COMPILER_WARNINGS=off \
    -D Kokkos_ENABLE_DEBUG=Off \
    -D Kokkos_ENABLE_AGGRESSIVE_VECTORIZATION=on \
    -D Kokkos_ENABLE_SERIAL=On \
    -D Kokkos_ENABLE_ATHREAD=on \
    -D Kokkos_ATHREAD_FUNCTOR_DIR="/home/export/online1/mdt00/shisuan/swiap/wjl/swKokkos-paper/swKokkos-AD/GCR/src" \
    -D MACHINE=${MACHINE}

    #-D CMAKE_FORTRAN_COMPILER="mpif90" \
    # -D Kokkos_ENABLE_ATHREAD=on \
    # -D Kokkos_ATHREAD_FUNCTOR_DIR="/home/export/base/shisuan/swiap/online/wjl/swKokkos-paper/dev/Kokkos-Athread/src"
fi

make -j
# make
