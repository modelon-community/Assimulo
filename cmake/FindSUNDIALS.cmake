# - Find SUNDIALS
# SUNDIALS, the SUite of Nonlinear and DIfferential/ALgebraic equation Solvers.
# https://computing.llnl.gov/projects/sundials
#
# The module defines the following variables:
#  SUNDIALS_VERSION, the version string
#  SUNDIALS_INCLUDE_DIRS, where to find sundials_dense.h, etc.
#  SUNDIALS_LIBRARIES, the libraries needed to use SUNDIALS
#  SUNDIALS_FOUND, If false, do not try to use SUNDIALS
# also defined, but not for general use are
#

find_path (SUNDIALS_INCLUDE_DIR sundials/sundials_config.h)

file (STRINGS ${SUNDIALS_INCLUDE_DIR}/sundials/sundials_config.h _VERSION_DEFINE_STRING REGEX "#define (SUNDIALS_VERSION|SUNDIALS_PACKAGE_VERSION) .*")
if (_VERSION_DEFINE_STRING)
  string (REGEX REPLACE "#define SUNDIALS[A-Z_]+VERSION \"([0-9\.]+)\"" "\\1" SUNDIALS_VERSION ${_VERSION_DEFINE_STRING})
endif ()

set(SUNDIALS_LIBRARIES)
set(SUNDIALS_COMPONENTS sunlinsoldense sunlinsolspgmr sunmatrixdense sunmatrixsparse core sunlinsolsuperlumt kinsol nvecserial cvode cvodes)
foreach (COMPONENT ${SUNDIALS_COMPONENTS})
  string(TOUPPER "${COMPONENT}" COMPONENT_UPPER)
  find_library (SUNDIALS_${COMPONENT_UPPER}_LIBRARY NAMES sundials_${COMPONENT})

  if (SUNDIALS_${COMPONENT_UPPER}_LIBRARY)
    list (APPEND SUNDIALS_LIBRARIES ${SUNDIALS_${COMPONENT_UPPER}_LIBRARY})
  endif ()
endforeach ()

set (SUNDIALS_INCLUDE_DIRS ${SUNDIALS_INCLUDE_DIR})

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(SUNDIALS REQUIRED_VARS SUNDIALS_INCLUDE_DIR SUNDIALS_CVODES_LIBRARY VERSION_VAR SUNDIALS_VERSION)

mark_as_advanced (
  SUNDIALS_LIBRARIES
  SUNDIALS_INCLUDE_DIR
  SUNDIALS_INCLUDE_DIRS)

if(NOT TARGET SUNDIALS::ALL)
  add_library(SUNDIALS::ALL UNKNOWN IMPORTED)
  set_target_properties(SUNDIALS::ALL
    PROPERTIES
    INTERFACE_INCLUDE_DIRECTORIES "${SUNDIALS_INCLUDE_DIRS}")
  target_link_libraries(SUNDIALS::ALL INTERFACE ${SUNDIALS_LIBRARIES})
endif()
