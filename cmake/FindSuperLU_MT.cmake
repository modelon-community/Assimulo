# - Find SuperLU_MT
# SuperLU_MT contains a set of subroutines to solve a sparse linear system
# https://github.com/xiaoyeli/superlu_mt
#
# The module defines the following variables:
#  SUPERLUMT_INCLUDE_DIRS, where to find superlu.h, etc.
#  SUPERLUMT_LIBRARIES, the libraries needed to use SuperLU_MT
#  SUPERLUMT_FOUND, If false, do not try to use SuperLU_MT
# also defined, but not for general use are
#  SUPERLUMT_LIBRARY, where to find the MPC library.
#

find_path (SUPERLUMT_INCLUDE_DIR slu_mt_cdefs.h
  PATH_SUFFIXES superlu_mt
)

find_library (SUPERLUMT_LIBRARY
  NAMES superlu_mt_OPENMP superlu_mt_THREAD
)

set (SUPERLUMT_LIBRARIES ${SUPERLUMT_LIBRARY})
set (SUPERLUMT_INCLUDE_DIRS ${SUPERLUMT_INCLUDE_DIR})

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(SuperLU_MT DEFAULT_MSG SUPERLUMT_LIBRARY SUPERLUMT_INCLUDE_DIRS)

mark_as_advanced (
  SUPERLUMT_LIBRARY
  SUPERLUMT_LIBRARIES
  SUPERLUMT_INCLUDE_DIR
  SUPERLUMT_INCLUDE_DIRS)



if(NOT TARGET SuperLU_MT::SuperLU_MT)
  add_library(SuperLU_MT::SuperLU_MT UNKNOWN IMPORTED)
  set_target_properties(SuperLU_MT::SuperLU_MT
    PROPERTIES
    IMPORTED_LOCATION ${SUPERLUMT_LIBRARY}
    INTERFACE_INCLUDE_DIRECTORIES "${SUPERLUMT_INCLUDE_DIR}")

  find_package(BLAS QUIET)
  if (BLAS_LIBRARIES)
    target_link_libraries(SuperLU_MT::SuperLU_MT INTERFACE ${BLAS_LIBRARIES})
  endif ()

  if (SUPERLUMT_LIBRARY MATCHES OPENMP)
    find_package(OpenMP QUIET)
    if (OpenMP_FOUND)
      target_link_libraries(SuperLU_MT::SuperLU_MT INTERFACE OpenMP::OpenMP_C)
    endif ()
    target_compile_definitions(SuperLU_MT::SuperLU_MT INTERFACE __OPENMP)
  endif ()

  find_library(MATH_LIBRARY NAMES m)
  if (MATH_LIBRARY)
    target_link_libraries(SuperLU_MT::SuperLU_MT INTERFACE ${MATH_LIBRARY})
  endif ()
endif()
