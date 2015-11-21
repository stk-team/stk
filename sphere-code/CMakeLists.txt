cmake_minimum_required(VERSION 2.8 FATAL_ERROR)
#set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake_modules/")
##############################################################################

project(sphereSamplingAnalysis)

set(LIBRARY_OUTPUT_PATH ./lib)
set(EXECUTABLE_OUTPUT_PATH .)

set(GARAGE_BIN_PATH bin)

file(MAKE_DIRECTORY ${GARAGE_BIN_PATH})

set(CMAKE_BUILD_TYPE Release)  # None | Debug  | Release | RelWithDebInfo | MinSizeRel

find_package(Boost 1.54.0 COMPONENTS system filesystem date_time program_options REQUIRED)
find_package(OpenMP)

if (OPENMP_FOUND)
        set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
        set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS} -std=gnu++0x -g -Wall -Wextra -Wno-dev -Wunused-variable")
else (OPENMP_FOUND)
        set(CMAKE_CXX_FLAGS "-std=gnu++0x -g -w -Wall -Wunused-variable")
endif(OPENMP_FOUND)

include_directories(
        ./src
        ${Boost_INCLUDE_DIR}
)

##sources
file(GLOB GARAGE_SOURCES
        /*.cpp
        io/*.cpp
        core/*.cpp
        sampler-sphere/*.cpp
        domaintransform/*.cpp
        sphere-analysis/*.cpp
        quaternions/*.cpp
        sampler-sphere/stk/*.cpp
)

##############################################################################
function(garage_make_exec exec_name)
        add_executable(${GARAGE_BIN_PATH}/${exec_name} ${GARAGE_SOURCES} ./demos/${exec_name}.cxx)
        target_link_libraries(
                ${GARAGE_BIN_PATH}/${exec_name}
                ${Boost_LIBRARIES}
        )
        install(PROGRAMS ${GARAGE_BIN_PATH}/${exec_name} DESTINATION bin)
endfunction()

##############################################################################

##Create executables by calling functions defined above:
#garage_make_exec(sampler-ccvt-sphere-demo1)
#garage_make_exec(variance-analysis-demo1-3)
#garage_make_exec(spherical-functions-demo1)
#garage_make_exec(samples-demo2)
garage_make_exec(powspec-sphere-demo2)
#garage_make_exec(powspec-projectivespace-hemisphere-to-sphere-demo1)
