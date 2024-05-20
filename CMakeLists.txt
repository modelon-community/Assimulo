cmake_minimum_required (VERSION 3.17)

project (assimulo LANGUAGES C)

list (APPEND CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake)

option(USE_FORTRAN "build Fortran modules" ON)
option(USE_SUNDIALS "build SUNDIALS modules" ON)
option(USE_SUPERLU "use SUPERLU lib" ON)
if (USE_FORTRAN)
  enable_language(Fortran)
endif ()

find_package (Python COMPONENTS Interpreter Development.Module NumPy REQUIRED)
include(GNUInstallDirs)
if (NOT DEFINED PYTHON_SITE_PACKAGES)
  if (WIN32)
    set (PYTHON_SITE_PACKAGES Lib/site-packages CACHE PATH "site-packages dir")
  else ()
    set (PYTHON_SITE_PACKAGES ${CMAKE_INSTALL_LIBDIR}/python${Python_VERSION_MAJOR}.${Python_VERSION_MINOR}/site-packages CACHE PATH "site-packages dir")
  endif ()
endif ()

find_package(LAPACK)

if (USE_SUPERLU)
  find_package (SuperLU_MT)
endif ()

if (USE_SUNDIALS)
  find_package (SUNDIALS CONFIG)

  if (SUNDIALS_FOUND)
    message(STATUS "Found SUNDIALS: ${SUNDIALS_DIR}  (found version \"${SUNDIALS_VERSION}\")")
  else ()
    # fallback to our module
    find_package (SUNDIALS MODULE)
    set(SUNDIALS_DIR ${SUNDIALS_INCLUDE_DIRS})
  endif ()

  if (SUNDIALS_FOUND)
    string (REGEX REPLACE "([0-9]+)\\..*" "\\1" SUNDIALS_MAJOR_VERSION "${SUNDIALS_VERSION}")
    string (REGEX REPLACE "[0-9]+\\.([0-9]+).*" "\\1" SUNDIALS_MINOR_VERSION "${SUNDIALS_VERSION}")
    string (REGEX REPLACE "[0-9]+\\.[0-9]+\\.([0-9]+).*" "\\1" SUNDIALS_PATCH_VERSION "${SUNDIALS_VERSION}")
    math(EXPR SUNDIALS_VERSION_NR "100000 * ${SUNDIALS_MAJOR_VERSION} + 100 * ${SUNDIALS_MINOR_VERSION} + ${SUNDIALS_PATCH_VERSION}")

    if (NOT DEFINED SUNDIALS_INCLUDE_DIR)
      find_path (SUNDIALS_INCLUDE_DIR sundials/sundials_config.h HINTS ${SUNDIALS_DIR}/../../../../include ${SUNDIALS_DIR}/../../../include)
    endif ()

    file (STRINGS ${SUNDIALS_INCLUDE_DIR}/sundials/sundials_config.h _SUNDIALS_INT64_DEFINE REGEX "^#define SUNDIALS_INT64_T")
    if (_SUNDIALS_INT64_DEFINE)
      set (SUNDIALS_VECTOR_SIZE 64)
    else ()
      set (SUNDIALS_VECTOR_SIZE 32)
    endif ()
    message(STATUS "SUNDIALS has vector type size ${SUNDIALS_VECTOR_SIZE} bits")

    # only SUNDIALS_SUPERLUMT_THREAD_TYPE is consistent across sundials versions
    file (STRINGS ${SUNDIALS_INCLUDE_DIR}/sundials/sundials_config.h _SUNDIALS_SUPERLUMT_DEFINE REGEX "^#define SUNDIALS_SUPERLUMT_THREAD_TYPE \"(PTHREAD|OPENMP)\"")
    if (_SUNDIALS_SUPERLUMT_DEFINE)
      set (SUNDIALS_WITH_SUPERLU True)
    else ()
      set (SUNDIALS_WITH_SUPERLU False)
    endif ()
    message(STATUS "SUNDIALS has SuperLU support: ${SUNDIALS_WITH_SUPERLU}")

    file (STRINGS ${SUNDIALS_INCLUDE_DIR}/sundials/sundials_config.h _SUNDIALS_CVODE_RTOL_VEC_DEFINE REGEX "^#define SUNDIALS_CVODE_RTOL_VEC")
    if (_SUNDIALS_CVODE_RTOL_VEC_DEFINE)
      set (SUNDIALS_CVODE_RTOL_VEC True)
    else ()
      set (SUNDIALS_CVODE_RTOL_VEC False)
    endif ()
    message(STATUS "SUNDIALS has CVode rtol vectors support: ${SUNDIALS_CVODE_RTOL_VEC}")

  endif ()
endif ()

set (cython_clones)
# copy /src tree into CMAKE_BINARY_DIR/assimulo
file (GLOB_RECURSE cython_sources src/*.pyx src/*.pxd src/*.py)
foreach (file ${cython_sources})
  if (EXISTS ${file})
    file (RELATIVE_PATH rel_file "${CMAKE_SOURCE_DIR}/src" "${file}")
    get_filename_component (rel_path ${rel_file} PATH)
    file (MAKE_DIRECTORY ${CMAKE_BINARY_DIR}/assimulo/${rel_path})
    set (cython_clone ${CMAKE_BINARY_DIR}/assimulo/${rel_file})
    execute_process(COMMAND ${CMAKE_COMMAND} -E copy ${file} ${cython_clone})
    add_custom_command (OUTPUT ${cython_clone}
      COMMENT "Copying ${file}"
      COMMAND ${CMAKE_COMMAND} -E copy ${file} ${cython_clone}
      DEPENDS ${file}
    )
    list (APPEND cython_clones ${cython_clone})
  endif ()
endforeach ()

# copy /thirdparty tree into CMAKE_BINARY_DIR
file (GLOB_RECURSE cython_sources thirdparty/*.pyx thirdparty/*.pxd thirdparty/*.py)
foreach (file ${cython_sources})
  if (EXISTS ${file})
    file (RELATIVE_PATH rel_file "${CMAKE_SOURCE_DIR}/thirdparty" "${file}")
    get_filename_component (rel_path ${rel_file} PATH)
    file (MAKE_DIRECTORY ${CMAKE_BINARY_DIR}/thirdparty/${rel_path})
    set (cython_clone ${CMAKE_BINARY_DIR}/thirdparty/${rel_file})
    execute_process(COMMAND ${CMAKE_COMMAND} -E copy ${file} ${cython_clone})
    add_custom_command (OUTPUT ${cython_clone}
      COMMENT "Copying ${file}"
      COMMAND ${CMAKE_COMMAND} -E copy ${file} ${cython_clone}
      DEPENDS ${file}
    )
    list (APPEND cython_clones ${cython_clone})
  endif ()
endforeach ()

macro(assimulo_add_cython_module pyx_file)
  cmake_parse_arguments(CMOD "" "DESTINATION" "SOURCES" ${ARGN})
  if (NOT DEFINED CMOD_DESTINATION)
    set (CMOD_DESTINATION assimulo)
  endif ()
  get_filename_component(name "${pyx_file}" NAME_WE)
  get_filename_component(subdir "${pyx_file}" DIRECTORY)

  # PYTHONNOUSERSITE=0 overrides skbuild and allows to use Cython from ~/.local when not available via system packages
  add_custom_command(
    OUTPUT ${name}.c
    COMMENT "Making ${name}.c from ${name}.pyx"
    COMMAND ${CMAKE_COMMAND} -E env PYTHONNOUSERSITE=0 ${Python_EXECUTABLE} -m cython -o ${name}.c
            --3str --fast-fail ${EXTRA_CYTHON_FLAGS}
            -I ${CMAKE_CURRENT_SOURCE_DIR}/src
            -I ${CMAKE_CURRENT_SOURCE_DIR}/src/lib
            ${CMAKE_CURRENT_BINARY_DIR}/${pyx_file}
    DEPENDS ${cython_clones} VERBATIM)

  python_add_library(${name} MODULE ${name}.c ${CMOD_SOURCES} WITH_SOABI)
  target_include_directories(${name} PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/src ${CMAKE_CURRENT_SOURCE_DIR}/${subdir})
  target_link_libraries(${name} PRIVATE Python::NumPy)
  target_compile_definitions(${name} PRIVATE NPY_NO_DEPRECATED_API=NPY_1_7_API_VERSION)
  if (CMAKE_C_COMPILER_ID MATCHES "GNU")
    target_compile_options(${name} PRIVATE $<$<COMPILE_LANGUAGE:C>:-O2 -fno-strict-aliasing>)
  endif ()
  install(TARGETS ${name} DESTINATION ${PYTHON_SITE_PACKAGES}/${CMOD_DESTINATION})
endmacro()

assimulo_add_cython_module(assimulo/explicit_ode.pyx
  SOURCES src/ode_event_locator.c)
assimulo_add_cython_module(assimulo/algebraic.pyx)
assimulo_add_cython_module(assimulo/implicit_ode.pyx)
assimulo_add_cython_module(assimulo/ode.pyx)
assimulo_add_cython_module(assimulo/problem.pyx)
assimulo_add_cython_module(assimulo/special_systems.pyx)
assimulo_add_cython_module(assimulo/support.pyx)
assimulo_add_cython_module(assimulo/solvers/euler.pyx
  DESTINATION assimulo/solvers)

if (SUNDIALS_FOUND)
  set(EXTRA_CYTHON_FLAGS -E SUNDIALS_VERSION_NR=${SUNDIALS_VERSION_NR} -E SUNDIALS_VECTOR_SIZE=${SUNDIALS_VECTOR_SIZE} -E SUNDIALS_WITH_SUPERLU=${SUNDIALS_WITH_SUPERLU} -E SUNDIALS_CVODE_RTOL_VEC=${SUNDIALS_CVODE_RTOL_VEC})

  assimulo_add_cython_module(assimulo/solvers/sundials.pyx
    DESTINATION assimulo/solvers)
  target_link_libraries(sundials PRIVATE SUNDIALS::cvodes_shared SUNDIALS::nvecserial_shared SUNDIALS::idas_shared)
  if (SUNDIALS_VERSION VERSION_GREATER_EQUAL 3)
    target_link_libraries(sundials PRIVATE SUNDIALS::sunlinsoldense_shared SUNDIALS::sunlinsolspgmr_shared SUNDIALS::sunmatrixdense SUNDIALS::sunmatrixsparse_shared)
  endif ()
  if (SUNDIALS_VERSION VERSION_GREATER_EQUAL 7)
    target_link_libraries(sundials PRIVATE SUNDIALS::core_shared)
  endif ()
  if (SUNDIALS_WITH_SUPERLU)
    target_link_libraries(sundials PRIVATE SuperLU_MT::SuperLU_MT)
    if (SUNDIALS_VERSION VERSION_GREATER_EQUAL 3)
      target_link_libraries(sundials PRIVATE SUNDIALS::sunlinsolsuperlumt_shared)
    endif ()
  endif ()

  assimulo_add_cython_module(assimulo/solvers/kinsol.pyx
    DESTINATION assimulo/solvers)
  target_link_libraries(kinsol PRIVATE SUNDIALS::kinsol SUNDIALS::nvecserial_shared)
  if (SUNDIALS_VERSION VERSION_GREATER_EQUAL 7)
    target_link_libraries(kinsol PRIVATE SUNDIALS::core_shared)
  endif ()
  if (SUNDIALS_WITH_SUPERLU)
    target_link_libraries(kinsol PRIVATE SuperLU_MT::SuperLU_MT)
  endif ()
endif ()


assimulo_add_cython_module(thirdparty/radau5/radau5ode.pyx
  SOURCES thirdparty/radau5/radau5.c thirdparty/radau5/radau5_io.c
  DESTINATION assimulo/lib)

if (SuperLU_MT_FOUND)
  target_sources(radau5ode PRIVATE thirdparty/radau5/superlu_double.c thirdparty/radau5/superlu_complex.c thirdparty/radau5/superlu_util.c)
  target_link_libraries(radau5ode PRIVATE SuperLU_MT::SuperLU_MT)
  target_compile_definitions(radau5ode PRIVATE __RADAU5_WITH_SUPERLU)
endif()

if (USE_FORTRAN)

  find_program(F2PY_EXECUTABLE NAMES f2py f2py3 REQUIRED)

  if (NOT DEFINED F2PY_INCLUDE_DIR)
    execute_process(COMMAND ${Python_EXECUTABLE} -c "from numpy import f2py; print(f2py.get_include())"
      OUTPUT_VARIABLE F2PY_INCLUDE_DIR
      OUTPUT_STRIP_TRAILING_WHITESPACE)
  endif()
  file(TO_CMAKE_PATH "${F2PY_INCLUDE_DIR}" F2PY_INCLUDE_DIR)
  message(STATUS "Using f2py include dir: ${F2PY_INCLUDE_DIR}")

  # object library with common f2py symbols
  add_library(fortranobject OBJECT ${F2PY_INCLUDE_DIR}/fortranobject.c)
  target_link_libraries(fortranobject PUBLIC Python::NumPy)
  target_compile_definitions(fortranobject PUBLIC NPY_NO_DEPRECATED_API=NPY_1_7_API_VERSION)
  target_include_directories(fortranobject PUBLIC ${F2PY_INCLUDE_DIR})
  set_target_properties(fortranobject PROPERTIES POSITION_INDEPENDENT_CODE TRUE)

  macro(assimulo_add_fortran_module name pyf_file)
    cmake_parse_arguments(FMOD "" "" "SOURCES" ${ARGN})
    add_custom_command(
      OUTPUT ${name}module.c
      COMMENT "Making ${name}module.c from ${pyf_file}"
      COMMAND ${F2PY_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/${pyf_file}
      DEPENDS ${pyf_file} VERBATIM)

    # may be created or not
    file(TOUCH ${CMAKE_CURRENT_BINARY_DIR}/${name}-f2pywrappers.f)
    file(TOUCH ${CMAKE_CURRENT_BINARY_DIR}/${name}-f2pywrappers2.f90)

    # separate lib for fortran
    add_library(${name}_fort STATIC
      ${FMOD_SOURCES}
      ${CMAKE_CURRENT_BINARY_DIR}/${name}-f2pywrappers.f
      ${CMAKE_CURRENT_BINARY_DIR}/${name}-f2pywrappers2.f90)
    set_target_properties(${name}_fort PROPERTIES POSITION_INDEPENDENT_CODE TRUE)
    set_source_files_properties(${CMAKE_CURRENT_BINARY_DIR}/${name}-f2pywrappers.f ${CMAKE_CURRENT_BINARY_DIR}/${name}-f2pywrappers2.f90 PROPERTIES GENERATED TRUE)
    if (CMAKE_Fortran_COMPILER_ID MATCHES "GNU")
      target_compile_options(${name}_fort PRIVATE $<$<COMPILE_LANGUAGE:Fortran>:-O2 --std=legacy>)
    endif ()

    python_add_library(${name} MODULE
      ${CMAKE_CURRENT_BINARY_DIR}/${name}module.c
      WITH_SOABI)
    set_source_files_properties(${CMAKE_CURRENT_BINARY_DIR}/${name}module.c PROPERTIES GENERATED TRUE)

    target_link_libraries(${name} PRIVATE fortranobject)
    target_link_libraries(${name} PRIVATE ${name}_fort)
    install(TARGETS ${name} DESTINATION ${PYTHON_SITE_PACKAGES}/assimulo/lib)
  endmacro()

  assimulo_add_fortran_module(dopri5 thirdparty/hairer/dopri5.pyf
    SOURCES thirdparty/hairer/dopri5.f)

  assimulo_add_fortran_module(rodas thirdparty/hairer/rodas_decsol.pyf
    SOURCES thirdparty/hairer/rodas_decsol.f)

  assimulo_add_fortran_module(radau5 thirdparty/hairer/radau_decsol.pyf
    SOURCES thirdparty/hairer/radau_decsol.f)

  assimulo_add_fortran_module(radar5 thirdparty/hairer/radar5.pyf
    SOURCES
      thirdparty/hairer/contr5.f90
      thirdparty/hairer/radar5_int.f90
      thirdparty/hairer/radar5.f90
      thirdparty/hairer/dontr5.f90
      thirdparty/hairer/decsol.f90
      thirdparty/hairer/dc_decdel.f90)

  assimulo_add_fortran_module(odepack thirdparty/odepack/odepack.pyf
    SOURCES
      thirdparty/odepack/opkdmain.f
      thirdparty/odepack/opkda1.f
      thirdparty/odepack/opkda2.f
      thirdparty/odepack/odepack_aux.f90)

  assimulo_add_fortran_module(odassl thirdparty/odassl/odassl.pyf
    SOURCES
      thirdparty/odassl/odassl.f
      thirdparty/odassl/odastp.f
      thirdparty/odassl/odacor.f
      thirdparty/odassl/odajac.f
      thirdparty/odassl/d1mach.f
      thirdparty/odassl/daxpy.f
      thirdparty/odassl/ddanrm.f
      thirdparty/odassl/ddatrp.f
      thirdparty/odassl/ddot.f
      thirdparty/odassl/ddwats.f
      thirdparty/odassl/dgefa.f
      thirdparty/odassl/dgesl.f
      thirdparty/odassl/dscal.f
      thirdparty/odassl/idamax.f
      thirdparty/odassl/xerrwv.f)
  if (CMAKE_Fortran_COMPILER_ID MATCHES "GNU")
    target_compile_options(odassl PRIVATE $<$<COMPILE_LANGUAGE:Fortran>:-fdefault-double-8 -fdefault-real-8>)
  endif ()

  assimulo_add_fortran_module(dasp3dp thirdparty/dasp3/dasp3dp.pyf
    SOURCES
      thirdparty/dasp3/DASP3.f
      thirdparty/dasp3/ANORM.f
      thirdparty/dasp3/CTRACT.f
      thirdparty/dasp3/DECOMP.f
      thirdparty/dasp3/HMAX.f
      thirdparty/dasp3/INIVAL.f
      thirdparty/dasp3/JACEST.f
      thirdparty/dasp3/PDERIV.f
      thirdparty/dasp3/PREPOL.f
      thirdparty/dasp3/SOLVE.f
      thirdparty/dasp3/SPAPAT.f)

  if (LAPACK_FOUND)
    assimulo_add_fortran_module(glimda thirdparty/glimda/glimda_complete.pyf
      SOURCES
        thirdparty/glimda/glimda_complete.f)
    target_link_libraries(glimda PRIVATE ${LAPACK_LIBRARIES})
  endif ()
endif ()

install(FILES src/__init__.py
  src/exception.py
  src/problem_algebraic.py
  src/algebraic.pxd
  src/explicit_ode.pxd
  src/implicit_ode.pxd
  src/ode.pxd
  src/problem.pxd
  src/support.pxd
  DESTINATION ${PYTHON_SITE_PACKAGES}/assimulo)
install(FILES src/lib/__init__.py
  src/lib/radau_core.py
  DESTINATION ${PYTHON_SITE_PACKAGES}/assimulo/lib)
install(FILES src/solvers/__init__.py
  src/solvers/dasp3.py
  src/solvers/glimda.py
  src/solvers/odassl.py
  src/solvers/odepack.py
  src/solvers/radar5.py
  src/solvers/radau5.py
  src/solvers/rosenbrock.py
  src/solvers/runge_kutta.py
  DESTINATION ${PYTHON_SITE_PACKAGES}/assimulo/solvers)
install(DIRECTORY examples
  DESTINATION ${PYTHON_SITE_PACKAGES}/assimulo)
