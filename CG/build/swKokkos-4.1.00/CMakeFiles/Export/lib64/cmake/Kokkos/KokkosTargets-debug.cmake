#----------------------------------------------------------------
# Generated CMake target import file for configuration "Debug".
#----------------------------------------------------------------

# Commands may need to know the format version.
set(CMAKE_IMPORT_FILE_VERSION 1)

# Import target "Kokkos::kokkoscore" for configuration "Debug"
set_property(TARGET Kokkos::kokkoscore APPEND PROPERTY IMPORTED_CONFIGURATIONS DEBUG)
set_target_properties(Kokkos::kokkoscore PROPERTIES
  IMPORTED_LINK_INTERFACE_LANGUAGES_DEBUG "C;CXX"
  IMPORTED_LOCATION_DEBUG "${_IMPORT_PREFIX}/lib64/libkokkoscore.a"
  )

list(APPEND _IMPORT_CHECK_TARGETS Kokkos::kokkoscore )
list(APPEND _IMPORT_CHECK_FILES_FOR_Kokkos::kokkoscore "${_IMPORT_PREFIX}/lib64/libkokkoscore.a" )

# Import target "Kokkos::kokkoscontainers" for configuration "Debug"
set_property(TARGET Kokkos::kokkoscontainers APPEND PROPERTY IMPORTED_CONFIGURATIONS DEBUG)
set_target_properties(Kokkos::kokkoscontainers PROPERTIES
  IMPORTED_LINK_INTERFACE_LANGUAGES_DEBUG "CXX"
  IMPORTED_LOCATION_DEBUG "${_IMPORT_PREFIX}/lib64/libkokkoscontainers.a"
  )

list(APPEND _IMPORT_CHECK_TARGETS Kokkos::kokkoscontainers )
list(APPEND _IMPORT_CHECK_FILES_FOR_Kokkos::kokkoscontainers "${_IMPORT_PREFIX}/lib64/libkokkoscontainers.a" )

# Import target "Kokkos::kokkossimd" for configuration "Debug"
set_property(TARGET Kokkos::kokkossimd APPEND PROPERTY IMPORTED_CONFIGURATIONS DEBUG)
set_target_properties(Kokkos::kokkossimd PROPERTIES
  IMPORTED_LINK_INTERFACE_LANGUAGES_DEBUG "CXX"
  IMPORTED_LOCATION_DEBUG "${_IMPORT_PREFIX}/lib64/libkokkossimd.a"
  )

list(APPEND _IMPORT_CHECK_TARGETS Kokkos::kokkossimd )
list(APPEND _IMPORT_CHECK_FILES_FOR_Kokkos::kokkossimd "${_IMPORT_PREFIX}/lib64/libkokkossimd.a" )

# Commands beyond this point should not need to know the version.
set(CMAKE_IMPORT_FILE_VERSION)
