########################################################################
project(stk)
########################################################################


cmake_minimum_required(VERSION 2.8)
cmake_policy(SET CMP0015 OLD)
SET(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${PROJECT_SOURCE_DIR}/modules")


set(CAIRO_ENABLED 1)
set(CMAKE_SHARED_MODULE_PREFIX)


message(STATUS "-------------------------------------------------------------------------------")
message(STATUS "STK Version 0.1")
message(STATUS "-------------------------------------------------------------------------------" )
message(STATUS "Source DIR is ${PROJECT_SOURCE_DIR}")
message(STATUS "Binary/Build DIR is ${PROJECT_BINARY_DIR}")
message(STATUS "Build type is ${CMAKE_BUILD_TYPE}")
message(STATUS "Installation prefix directory is " ${CMAKE_INSTALL_PREFIX})
message(STATUS "Host system is " ${CMAKE_HOST_SYSTEM} " with processor " ${CMAKE_HOST_SYSTEM_PROCESSOR})
message(STATUS "Target system is " ${CMAKE_SYSTEM} " with processor " ${CMAKE_SYSTEM_PROCESSOR})
message(STATUS "-------------------------------------------------------------------------------" )
message(STATUS "------------------------------------DEPS---------------------------------------" )

#------------------------------------------------------------------------------
# Offer the user the choice of overriding the installation directories
#------------------------------------------------------------------------------
set(INSTALL_LIB_DIR lib CACHE PATH "Installation directory for libraries")
set(INSTALL_BIN_DIR bin CACHE PATH "Installation directory for executables")
set(INSTALL_INCLUDE_DIR include CACHE PATH "Installation directory for header files")
set(INSTALL_DATA_DIR share CACHE PATH "Installation directory for data files")

########################################################################
## DEPS
########################################################################

find_package(PNG REQUIRED)
find_package(Boost 1.42.0 COMPONENTS thread program_options unit_test_framework timer REQUIRED)
find_package(FFTW REQUIRED)
find_package(CGAL REQUIRED COMPONENTS Core)
find_package(GMP REQUIRED)

if(CAIRO_ENABLED)
  find_package(Cairo REQUIRED)
  add_definitions(-DCAIRO_ENABLED)
endif(CAIRO_ENABLED)

add_definitions(-DDIR_MODULE="${STK_MODULE_PATH}/")

## disabling for now.. not supported on clang
#set(CMAKE_CXX_FLAGS ${CMAKE_CXX_FLAGS} -frounding-math)

include_directories(
  ${PROJECT_SOURCE_DIR}/src
  ${PNG_INCLUDE_DIR}
  ${Boost_INCLUDE_DIR}
  ${FFTW_INCLUDE_DIR}
  ${GMP_INCLUDE_DIR}
  ${CGAL_INCLUDE_DIRS}
  ${CGAL_3RD_PARTY_INCLUDE_DIRS}
  /opt/local/include/cairo/
  /usr/include/suitesparse
  )

link_directories(
  ${PROJECT_SOURCE_DIR}/lib
  ${Boost_LIBRARY_DIRS}
  ${FFTW_INCLUDE_DIR}
  ${CGAL_LIBRARIES_DIR}
  ${CGAL_3RD_PARTY_LIBRARIES_DIRS}
  ${CAIRO_LIBRARY_DIRS}
  )

set(STK_DEPS_LIBRARIES
  ${CGAL_LIBRARIES}
  ${CGAL_CORE_LIBRARY}
  ${CGAL_3RD_PARTY_LIBRARIES}
  ${PNG_LIBRARIES}
  ${FFTW_LIBRARIES}
  ${CAIRO_LIBRARIES}
  ${GMPXX_LIBRARIES}
  ${GMP_LIBRARIES}
  ${Boost_LIBRARIES})

include(${CGAL_USE_FILE})

message(STATUS "-------------------------------------------------------------------------------" )

########################################################################
# SOURCES
########################################################################

file(GLOB stk_sources
  src/stk/*.cpp
  src/stk/conversion/*.cpp
  src/stk/diff-space/*.cpp
  src/stk/fourier/*.cpp
  src/stk/function/*.cpp
  src/stk/io/*.cpp
  src/stk/optimization/*.cpp
  src/stk/sampler/*.cpp
  src/stk/unstable/*.cpp
  )

file(GLOB stk_headers src/stk/*.hpp)
file(GLOB stk_headers_conversion src/stk/conversion/*.hpp)
file(GLOB stk_headers_diff-space src/stk/diff-space/*.hpp)
file(GLOB stk_headers_fourier src/stk/fourier/*.hpp)
file(GLOB stk_headers_function src/stk/function/*.hpp)
file(GLOB stk_headers_io src/stk/io/*.hpp)
file(GLOB stk_headers_optimization src/stk/optimization/*.hpp)
file(GLOB stk_headers_sampler src/stk/sampler/*.hpp)
file(GLOB stk_headers_unstable src/stk/unstable/*.hpp)

########################################################################
# LIBRARIES
########################################################################

add_library(stk SHARED ${stk_sources})
add_library(stk_static STATIC ${stk_sources})

target_link_libraries(
  stk
  ${STK_DEPS_LIBRARIES}
  )
target_link_libraries(
  stk_static
  ${STK_DEPS_LIBRARIES}
  )

set_target_properties(stk PROPERTIES VERSION 0.1)
set_target_properties(stk PROPERTIES LINKER_LANGUAGE CXX)
set_target_properties(stk_static PROPERTIES VERSION 0.1)
set_target_properties(stk_static PROPERTIES LINKER_LANGUAGE CXX)

install(TARGETS stk stk_static
  DESTINATION lib
  )

install(FILES ${stk_headers}
  DESTINATION ${INSTALL_INCLUDE_DIR}/stk)
install(FILES ${stk_headers_conversion}
  DESTINATION ${INSTALL_INCLUDE_DIR}/stk/conversion)
install(FILES ${stk_headers_fourier}
  DESTINATION ${INSTALL_INCLUDE_DIR}/stk/fourier)
install(FILES ${stk_headers_function}
  DESTINATION ${INSTALL_INCLUDE_DIR}/stk/function)
install(FILES ${stk_headers_io}
  DESTINATION ${INSTALL_INCLUDE_DIR}/stk/io)
install(FILES ${stk_headers_optimization}
  DESTINATION ${INSTALL_INCLUDE_DIR}/stk/optimization)
install(FILES ${stk_headers_sampler}
  DESTINATION ${INSTALL_INCLUDE_DIR}/stk/sampler)
install(FILES ${stk_headers_unstable}
  DESTINATION ${INSTALL_INCLUDE_DIR}/stk/unstable)

install(FILES data/cmake/FindSTK.cmake
  DESTINATION /usr/share/cmake-2.8/Modules)

if(CUDA_ENABLED)
  install(FILES data/cmake/enabled/STKCuda.cmake
    DESTINATION /usr/local/share/stk/cmake)
else(CUDA_ENABLED)
  install(FILES data/cmake/disabled/STKCuda.cmake
    DESTINATION /usr/local/share/stk/cmake)
endif(CUDA_ENABLED)

if(CAIRO_ENABLED)
  install(FILES data/cmake/enabled/STKCairo.cmake
    DESTINATION /usr/local/share/stk/cmake)
else(CAIRO_ENABLED)
  install(FILES data/cmake/disabled/STKCairo.cmake
    DESTINATION /usr/local/share/stk/cmake)
endif(CAIRO_ENABLED)

########################################################################
# EXECUTABLES
########################################################################

function(stk_make_program prog_name prog_libs)
  add_executable(stk-${prog_name} src/${prog_name}.cpp)
  target_link_libraries(
    stk-${prog_name}
    stk
    ${STK_DEPS_LIBRARIES}
    ${prog_libs}
    )
  install(PROGRAMS stk-${prog_name} DESTINATION bin)
endfunction()

stk_make_program(fourier "")
stk_make_program(integration-test "")
stk_make_program(pts-info "")
stk_make_program(discrepancy "")
#~ stk_make_program(checkerboard "")
stk_make_program(zoneplate "")
#~ stk_make_program(zoneplate2 "")
#~ stk_make_program(zoneplate-variance -ldl)
#~ stk_make_program(lloyd "")
#~ stk_make_program(convert "")
#~ stk_make_program(bpts "")
#~ stk_make_program(ddf "")
#~ stk_make_program(integration -ldl)

#~ stk_make_program(discrepancy3d "")
#~ stk_make_program(voronoi "")
#~ stk_make_program(capacity "")
#~ stk_make_program(cdf "")
#~ stk_make_program(discrepancy4d "")

#~ stk_make_program(pts-merge "")
#~ stk_make_program(pts-take "")

if(CAIRO_ENABLED)
  stk_make_program(draw-pts "")
  #~ stk_make_program(draw-histo "")
endif(CAIRO_ENABLED)

########################################################################
# SAMPLER
########################################################################

stk_make_program(sampler-whitenoiseT2 "")
stk_make_program(sampler-stratifiedT2 "")
stk_make_program(sampler-sobolT2 "")
stk_make_program(sampler-haltonT2 "")

add_executable(stk-sampler-poissondiskT2 src/sampler-poissondisk/sampler-poissondiskT2.cpp)
target_link_libraries(
  stk-sampler-poissondiskT2
  stk
  ${STK_DEPS_LIBRARIES}
  ${prog_libs}
  )
install(PROGRAMS stk-sampler-poissondiskT2 DESTINATION bin)

add_executable(stk-sampler-fpoT2 src/sampler-fpo/sampler-fpoT2.cpp src/sampler-fpo/DT.cpp)
target_link_libraries(
  stk-sampler-fpoT2
  stk
  ${STK_DEPS_LIBRARIES}
  ${prog_libs}
  )
install(PROGRAMS stk-sampler-poissondiskT2 DESTINATION bin)

#~ add_library(sampler/bnot MODULE
#~ add_executable(${STK_BIN_PATH}/bnot
#~ src/sampler/bnot/bnot.cpp
#~ src/sampler/bnot/matrix/sparse_array.cpp
#~ src/sampler/bnot/matrix/sparse_matrix.cpp
#~ src/sampler/bnot/energy.cpp
#~ src/sampler/bnot/init.cpp
#~ src/sampler/bnot/optimizer.cpp
#~ src/sampler/bnot/regularity.cpp
#~ src/sampler/bnot/sites.cpp)

#~ target_link_libraries(
#~ sampler/bnot
#~ ${CGAL_LIBRARIES} -lgmp
#~ -lCGAL -lCGAL_Core 
#~ stk
#~ -lm -lamd -lcamd -lcolamd -lccolamd -lcholmod -lspqr -ltbb -lmetis -lblas -llapack
#~ )
#~ install(TARGETS sampler/bnot DESTINATION ${STK_MODULE_PATH}/samplers)


########################################################################
# UNIT TEST
########################################################################

enable_testing()

function(stk_make_test test_name)
  add_executable(${STK_TEST_BIN_PATH}/test-${test_name} unit-test/src/${test_name}.cpp)
  target_link_libraries(
    ${STK_TEST_BIN_PATH}/test-${test_name}
    stk
    ${STK_DEPS_LIBRARIES}
    ${Boost_LIBRARIES}
    ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY_DEBUG}
    )
  add_test(test-${test_name} ${STK_TEST_BIN_PATH}/test-${test_name})
endfunction()

#~ stk_make_test(vector)
#~ stk_make_test(domain)
#~ stk_make_test(pointset)
#~ stk_make_test(pointset-stream)
#~ stk_make_test(array)

