IF(KOKKOS_ENABLE_ATHREAD)
KOKKOS_INCLUDE_DIRECTORIES(
  ${CMAKE_CURRENT_BINARY_DIR}
  ${CMAKE_CURRENT_SOURCE_DIR}
  "${CMAKE_CURRENT_SOURCE_DIR}/Athread"
  ${KOKKOS_TOP_BUILD_DIR}
)
ELSE()
KOKKOS_INCLUDE_DIRECTORIES(
  ${CMAKE_CURRENT_BINARY_DIR}
  ${CMAKE_CURRENT_SOURCE_DIR}
  ${KOKKOS_TOP_BUILD_DIR}
)
ENDIF()

IF (NOT desul_FOUND)
  IF(KOKKOS_ENABLE_CUDA)
    SET(DESUL_ATOMICS_ENABLE_CUDA ON)
  ENDIF()
  IF(KOKKOS_ENABLE_CUDA_RELOCATABLE_DEVICE_CODE)
    SET(DESUL_ATOMICS_ENABLE_CUDA_SEPARABLE_COMPILATION ON)
  ENDIF()
  IF(KOKKOS_ENABLE_HIP)
    SET(DESUL_ATOMICS_ENABLE_HIP ON)
  ENDIF()
  IF(KOKKOS_ENABLE_HIP_RELOCATABLE_DEVICE_CODE)
    SET(DESUL_ATOMICS_ENABLE_HIP_SEPARABLE_COMPILATION ON)
  ENDIF()
  IF(KOKKOS_ENABLE_SYCL)
    SET(DESUL_ATOMICS_ENABLE_SYCL ON)
  ENDIF()
  IF(KOKKOS_ENABLE_OPENMPTARGET)
    SET(DESUL_ATOMICS_ENABLE_OPENMP ON) # not a typo Kokkos OpenMPTarget -> Desul OpenMP
  ENDIF()
  IF(KOKKOS_ENABLE_ATHREAD)
    # TODO
    SET(DESUL_ATOMICS_ENABLE_ATHREAD OFF)
  ENDIF()
  CONFIGURE_FILE(
    ${CMAKE_CURRENT_SOURCE_DIR}/../../tpls/desul/Config.hpp.cmake.in
    ${CMAKE_CURRENT_BINARY_DIR}/desul/atomics/Config.hpp
  )
  KOKKOS_INCLUDE_DIRECTORIES(
    ${CMAKE_CURRENT_SOURCE_DIR}/../../tpls/desul/include
  )
ENDIF()

INSTALL (DIRECTORY
  "${CMAKE_CURRENT_SOURCE_DIR}/"
  DESTINATION ${KOKKOS_HEADER_DIR}
  FILES_MATCHING
  PATTERN "*.hpp"
  PATTERN "*.h"
)

SET(KOKKOS_CORE_SRCS)
APPEND_GLOB(KOKKOS_CORE_SRCS ${CMAKE_CURRENT_SOURCE_DIR}/impl/*.cpp)
SET(KOKKOS_CORE_HEADERS)
APPEND_GLOB(KOKKOS_CORE_HEADERS ${CMAKE_CURRENT_SOURCE_DIR}/*.hpp)
APPEND_GLOB(KOKKOS_CORE_HEADERS ${CMAKE_CURRENT_SOURCE_DIR}/impl/*.hpp)

IF (KOKKOS_ENABLE_CUDA)
  APPEND_GLOB(KOKKOS_CORE_SRCS ${CMAKE_CURRENT_SOURCE_DIR}/Cuda/*.cpp)
  APPEND_GLOB(KOKKOS_CORE_HEADERS ${CMAKE_CURRENT_SOURCE_DIR}/Cuda/*.hpp)
ENDIF()

IF (KOKKOS_ENABLE_OPENMP)
  APPEND_GLOB(KOKKOS_CORE_SRCS ${CMAKE_CURRENT_SOURCE_DIR}/OpenMP/*.cpp)
  APPEND_GLOB(KOKKOS_CORE_HEADERS ${CMAKE_CURRENT_SOURCE_DIR}/OpenMP/*.hpp)
ENDIF()

IF (KOKKOS_ENABLE_OPENMPTARGET)
  APPEND_GLOB(KOKKOS_CORE_SRCS ${CMAKE_CURRENT_SOURCE_DIR}/OpenMPTarget/*.cpp)
  APPEND_GLOB(KOKKOS_CORE_HEADERS ${CMAKE_CURRENT_SOURCE_DIR}/OpenMPTarget/*.hpp)
ENDIF()

IF (KOKKOS_ENABLE_OPENACC)
  APPEND_GLOB(KOKKOS_CORE_SRCS ${CMAKE_CURRENT_SOURCE_DIR}/OpenACC/*.cpp)
  APPEND_GLOB(KOKKOS_CORE_HEADERS ${CMAKE_CURRENT_SOURCE_DIR}/OpenACC/*.hpp)
ENDIF()

IF (KOKKOS_ENABLE_THREADS)
  APPEND_GLOB(KOKKOS_CORE_SRCS ${CMAKE_CURRENT_SOURCE_DIR}/Threads/*.cpp)
  APPEND_GLOB(KOKKOS_CORE_HEADERS ${CMAKE_CURRENT_SOURCE_DIR}/Threads/*.hpp)
ENDIF()

IF (KOKKOS_ENABLE_ATHREAD)
  execute_process(COMMAND python 
    ${CMAKE_CURRENT_SOURCE_DIR}/Athread/Kokkos_Athread_Register_Template.py 
    ${Kokkos_ATHREAD_FUNCTOR_DIR}
    ${CMAKE_CURRENT_SOURCE_DIR}/Athread/Kokkos_Athread_KernelLaunch_Slave.cpp)

  APPEND_GLOB(KOKKOS_CORE_HEADERS ${CMAKE_CURRENT_SOURCE_DIR}/Athread/*.h)
  APPEND_GLOB(KOKKOS_CORE_HEADERS ${CMAKE_CURRENT_SOURCE_DIR}/Athread/*.hpp)


  APPEND_GLOB(KOKKOS_CORE_SRCS 
      ${CMAKE_CURRENT_SOURCE_DIR}/Athread/Kokkos_Athread_HashTable_Host.cpp)
  APPEND_GLOB(KOKKOS_CORE_SRCS 
      ${CMAKE_CURRENT_SOURCE_DIR}/Athread/Kokkos_Athread_Task.cpp)
  APPEND_GLOB(KOKKOS_CORE_SRCS 
      ${CMAKE_CURRENT_SOURCE_DIR}/Athread/Kokkos_Athread.cpp)

  APPEND_GLOB(KOKKOS_CORE_SRCS_C 
      ${CMAKE_CURRENT_SOURCE_DIR}/Athread/Kokkos_Athread_KernelLaunch.c)

  APPEND_GLOB(KOKKOS_CORE_SRCS_C 
	  ${CMAKE_CURRENT_SOURCE_DIR}/Athread/Kokkos_Athread_FastThreadSpawn.c)
  APPEND_GLOB(KOKKOS_CORE_SRCS ${KOKKOS_CORE_SRCS_C})

  set_source_files_properties(${KOKKOS_CORE_SRCS} PROPERTIES COMPILE_FLAGS "-mhost -fPIC -msimd -faddress_align=64")

  # /usr/sw/yyzlib/xMath-SACA
  APPEND_GLOB(KOKKOS_CORE_SRCS_SLAVE 
      ${CMAKE_CURRENT_SOURCE_DIR}/Athread/Kokkos_Athread_KernelLaunch_Slave.cpp)
  APPEND_GLOB(KOKKOS_CORE_SRCS_SLAVE 
      ${CMAKE_CURRENT_SOURCE_DIR}/Athread/Kokkos_Athread_Utils_Slave.cpp)
  APPEND_GLOB(KOKKOS_CORE_SRCS_SLAVE 
      ${CMAKE_CURRENT_SOURCE_DIR}/Athread/Kokkos_Athread_HashTable_Slave.cpp)
  APPEND_GLOB(KOKKOS_CORE_SRCS_SLAVE 
      ${CMAKE_CURRENT_SOURCE_DIR}/Athread/Kokkos_Athread_FastThreadSpawn_Slave.c)
  set_source_files_properties(${KOKKOS_CORE_SRCS_SLAVE} 
      PROPERTIES COMPILE_FLAGS "-mslave -faddress_align=64 -fPIC -msimd")

  #APPEND_GLOB(KOKKOS_CORE_SRCS_SLAVE_C 
  #    ${CMAKE_CURRENT_SOURCE_DIR}/Athread/Kokkos_Athread_MathematicalFunctionsSlave.c)
  APPEND_GLOB(KOKKOS_CORE_SRCS_SLAVE_C 
      ${CMAKE_CURRENT_SOURCE_DIR}/Athread/Kokkos_Athread_KernelLaunch_SlaveC.c)
  set_source_files_properties(${KOKKOS_CORE_SRCS_SLAVE_C} 
      PROPERTIES COMPILE_FLAGS "-mftz -mieee -mslave -msimd -faddress_align=64 -fmath-opt=sin,exp,cos,sqrt -lm_slave -fPIC -O3")

  APPEND_GLOB(KOKKOS_CORE_SRCS ${KOKKOS_CORE_SRCS_SLAVE})
  APPEND_GLOB(KOKKOS_CORE_SRCS ${KOKKOS_CORE_SRCS_SLAVE_C})
  
ENDIF()

IF (KOKKOS_ENABLE_HIP)
  APPEND_GLOB(KOKKOS_CORE_SRCS ${CMAKE_CURRENT_SOURCE_DIR}/HIP/*.cpp)
  APPEND_GLOB(KOKKOS_CORE_HEADERS ${CMAKE_CURRENT_SOURCE_DIR}/HIP/*.hpp)
ENDIF()

IF (KOKKOS_ENABLE_HPX)
  APPEND_GLOB(KOKKOS_CORE_SRCS ${CMAKE_CURRENT_SOURCE_DIR}/HPX/*.cpp)
  APPEND_GLOB(KOKKOS_CORE_HEADERS ${CMAKE_CURRENT_SOURCE_DIR}/HPX/*.hpp)
ENDIF()

IF (NOT KOKKOS_ENABLE_MEMKIND)
  LIST(REMOVE_ITEM KOKKOS_CORE_SRCS ${CMAKE_CURRENT_SOURCE_DIR}/impl/Kokkos_HBWSpace.cpp)
ENDIF()

IF (KOKKOS_ENABLE_SERIAL)
  APPEND_GLOB(KOKKOS_CORE_SRCS ${CMAKE_CURRENT_SOURCE_DIR}/Serial/*.cpp)
  APPEND_GLOB(KOKKOS_CORE_HEADERS ${CMAKE_CURRENT_SOURCE_DIR}/Serial/*.hpp)
ENDIF()

IF (KOKKOS_ENABLE_SYCL)
  APPEND_GLOB(KOKKOS_CORE_SRCS ${CMAKE_CURRENT_SOURCE_DIR}/SYCL/*.cpp)
  APPEND_GLOB(KOKKOS_CORE_HEADERS ${CMAKE_CURRENT_SOURCE_DIR}/SYCL/*.hpp)
ENDIF()

IF (NOT desul_FOUND)
  IF (KOKKOS_ENABLE_CUDA)
    APPEND_GLOB(KOKKOS_CORE_SRCS ${CMAKE_CURRENT_SOURCE_DIR}/../../tpls/desul/src/Lock_Array_CUDA.cpp)
  ELSEIF (KOKKOS_ENABLE_HIP)
    APPEND_GLOB(KOKKOS_CORE_SRCS ${CMAKE_CURRENT_SOURCE_DIR}/../../tpls/desul/src/Lock_Array_HIP.cpp)
  ELSEIF (KOKKOS_ENABLE_SYCL)
    APPEND_GLOB(KOKKOS_CORE_SRCS ${CMAKE_CURRENT_SOURCE_DIR}/../../tpls/desul/src/Lock_Array_SYCL.cpp)
  ENDIF()
  APPEND_GLOB(KOKKOS_CORE_HEADERS ${CMAKE_CURRENT_SOURCE_DIR}/../../tpls/desul/include/desul/*.hpp)
  APPEND_GLOB(KOKKOS_CORE_HEADERS ${CMAKE_CURRENT_SOURCE_DIR}/../../tpls/desul/include/desul/*/*.hpp)
  APPEND_GLOB(KOKKOS_CORE_HEADERS ${CMAKE_CURRENT_SOURCE_DIR}/../../tpls/desul/include/desul/*/*/*.hpp)
  APPEND_GLOB(KOKKOS_CORE_HEADERS ${CMAKE_CURRENT_SOURCE_DIR}/../../tpls/desul/include/*/*/*.inc*)
  APPEND_GLOB(KOKKOS_CORE_HEADERS ${CMAKE_CURRENT_BINARY_DIR}/desul/*.hpp)

  INSTALL (DIRECTORY
    "${CMAKE_CURRENT_SOURCE_DIR}/../../tpls/desul/include/desul"
    "${CMAKE_CURRENT_BINARY_DIR}/desul"
    DESTINATION ${KOKKOS_HEADER_DIR}
    FILES_MATCHING
    PATTERN "*.inc"
    PATTERN "*.inc_*"
    PATTERN "*.hpp"
  )

  MESSAGE(STATUS "Using internal desul_atomics copy")
ELSE()
  MESSAGE(STATUS "Using external desul_atomics install found at:")
  MESSAGE(STATUS "  " ${desul_DIR})
ENDIF()


KOKKOS_ADD_LIBRARY(
  kokkoscore
  SOURCES ${KOKKOS_CORE_SRCS}
  HEADERS ${KOKKOS_CORE_HEADERS}
  ADD_BUILD_OPTIONS # core should be given all the necessary compiler/linker flags
)
IF (KOKKOS_ENABLE_ATHREAD)
  set_source_files_properties(${KOKKOS_CORE_SRCS_C} PROPERTIES LANGUAGE C)
  set_source_files_properties(${KOKKOS_CORE_SRCS_SLAVE_C} PROPERTIES LANGUAGE C)
ENDIF()

KOKKOS_LIB_INCLUDE_DIRECTORIES(kokkoscore
  ${KOKKOS_TOP_BUILD_DIR}
  ${CMAKE_CURRENT_BINARY_DIR}
  ${CMAKE_CURRENT_SOURCE_DIR}
)
IF (NOT desul_FOUND)
  KOKKOS_LIB_INCLUDE_DIRECTORIES(kokkoscore
    ${CMAKE_CURRENT_SOURCE_DIR}/../../tpls/desul/include
  )
ENDIF()

IF (Kokkos_ENABLE_IMPL_MDSPAN)
  MESSAGE(STATUS "Experimental mdspan support is enabled")

  # Some compilers now include mdspan... we just flag on their version
  # for now until we can get some compiler detection support
  include(CheckIncludeFileCXX)
  check_include_file_cxx(experimental/mdspan KOKKOS_COMPILER_SUPPORTS_EXPERIMENTAL_MDSPAN)
  check_include_file_cxx(mdspan KOKKOS_COMPILER_SUPPORTS_MDSPAN)

  if (Kokkos_ENABLE_MDSPAN_EXTERNAL)
    MESSAGE(STATUS "Using external mdspan")
    target_link_libraries(kokkoscore PUBLIC std::mdspan)
  elseif(KOKKOS_COMPILER_SUPPORTS_MDSPAN AND NOT Kokkos_ENABLE_IMPL_SKIP_COMPILER_MDSPAN)
    message(STATUS "Using compiler-supplied mdspan")
  elseif(KOKKOS_COMPILER_SUPPORTS_EXPERIMENTAL_MDSPAN AND NOT Kokkos_ENABLE_IMPL_SKIP_COMPILER_MDSPAN)
    message(STATUS "Using compiler-supplied experimental/mdspan")
  else()
    KOKKOS_LIB_INCLUDE_DIRECTORIES(
      kokkoscore
      ${CMAKE_CURRENT_SOURCE_DIR}/../../tpls/mdspan/include
    )

    APPEND_GLOB(KOKKOS_CORE_HEADERS ${CMAKE_CURRENT_SOURCE_DIR}/../../tpls/mdspan/include/experimental/__p0009_bits/*.hpp)
    APPEND_GLOB(KOKKOS_CORE_HEADERS ${CMAKE_CURRENT_SOURCE_DIR}/../../tpls/mdspan/include/experimental/mdspan)

    INSTALL (DIRECTORY
      "${CMAKE_CURRENT_SOURCE_DIR}/../../tpls/mdspan/include/"
      DESTINATION ${KOKKOS_HEADER_DIR}
      FILES_MATCHING
      PATTERN "mdspan"
      PATTERN "*.hpp"
    )
    MESSAGE(STATUS "Using internal mdspan directory ${CMAKE_CURRENT_SOURCE_DIR}/../../tpls/mdspan/include")
  endif()
ENDIF()

KOKKOS_LINK_TPL(kokkoscore PUBLIC HWLOC)
KOKKOS_LINK_TPL(kokkoscore PUBLIC MEMKIND)
IF (NOT KOKKOS_ENABLE_COMPILE_AS_CMAKE_LANGUAGE)
  KOKKOS_LINK_TPL(kokkoscore PUBLIC CUDA)
ENDIF()
KOKKOS_LINK_TPL(kokkoscore PUBLIC HPX)
KOKKOS_LINK_TPL(kokkoscore PUBLIC LIBDL)
KOKKOS_LINK_TPL(kokkoscore PUBLIC LIBRT)
# On *nix-like systems (Linux, macOS) we need pthread for C++ std::thread
IF (NOT WIN32)
  KOKKOS_LINK_TPL(kokkoscore PUBLIC THREADS)
ENDIF()
IF (NOT KOKKOS_ENABLE_COMPILE_AS_CMAKE_LANGUAGE)
  KOKKOS_LINK_TPL(kokkoscore PUBLIC ROCM)
  KOKKOS_LINK_TPL(kokkoscore PUBLIC ONEDPL)
ENDIF()

# FIXME: We need a proper solution to figure out whether to enable
#        libatomic
# Most compilers only require libatomic for 128-bit CAS
# I (CT) had removed 128bit CAS from desul to not need libatomic.
IF (KOKKOS_ENABLE_OPENMPTARGET)
  target_link_libraries(kokkoscore PUBLIC atomic)
ENDIF()

IF (desul_FOUND)
  target_link_libraries(kokkoscore PUBLIC desul_atomics)
ENDIF()

# FIXME_TRILINOS Trilinos doesn't allow for Kokkos to use find_dependency so we
# just append the flags in cmake/kokkos_tpls.cmake instead of linking with the
# OpenMP target.
IF(Kokkos_ENABLE_OPENMP AND NOT KOKKOS_HAS_TRILINOS)
  target_link_libraries(kokkoscore PUBLIC OpenMP::OpenMP_CXX)
ENDIF()

IF(Kokkos_ENABLE_ATHREAD)
  target_link_libraries(kokkoscore PUBLIC -lstdc++ -lm_slave -lm -Wl,--as-needed -mhybrid -static)
  # target_link_libraries(kokkoscore PUBLIC -mhybrid -static)
  # target_link_libraries(kokkoscore PUBLIC -mdynamic -shared)
ENDIF()

KOKKOS_LINK_TPL(kokkoscore PUBLIC LIBQUADMATH)
