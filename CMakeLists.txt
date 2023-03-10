cmake_minimum_required(VERSION 3.3)

project(SequentialLineSearch CXX)

if(SEQUENTIAL_LINE_SEARCH_BUILD_PYTHON_BINDING)
	set(CMAKE_POSITION_INDEPENDENT_CODE ON)
endif()

set(CMAKE_CXX_STANDARD 11)

# Append possible paths to Qt5 for macOS
if(APPLE AND EXISTS /usr/local/opt/qt5)
	list(APPEND CMAKE_PREFIX_PATH "/usr/local/opt/qt5")
endif()
if(APPLE AND EXISTS /usr/local/opt/qt)
	list(APPEND CMAKE_PREFIX_PATH "/usr/local/opt/qt")
endif()
if(APPLE AND EXISTS /opt/homebrew/opt/qt5)
	list(APPEND CMAKE_PREFIX_PATH "/opt/homebrew/opt/qt5")
endif()

option(SEQUENTIAL_LINE_SEARCH_BUILD_COMMAND_DEMOS                 "" ON )
option(SEQUENTIAL_LINE_SEARCH_BUILD_VISUAL_DEMOS                  "" OFF)
option(SEQUENTIAL_LINE_SEARCH_BUILD_PHOTO_DEMOS                   "" OFF)
option(SEQUENTIAL_LINE_SEARCH_BUILD_PYTHON_BINDING                "" OFF)
option(SEQUENTIAL_LINE_SEARCH_USE_NOISELESS_FORMULATION           "" OFF)
option(SEQUENTIAL_LINE_SEARCH_USE_PARALLELIZED_MULTI_START_SEARCH "" OFF)

###############################################
# Handling external libraries
###############################################

# Eigen
IF (WIN32)
	# message("${EIGEN3_INCLUDE_DIR}")
	# You have to set these variables as Windows environment variables:
    # EIGEN3_INCLUDE_DIR <- %EIGEN3_ROOT%
    # EIGEN3_DIR <- %EIGEN3_ROOT%\cmake
    #
    # EIGEN3_INCLUDE_DIR: variable needed for file %EIGEN3_ROOT%/cmake/FindEigen3.cmake
    #
    # CMAKE_MODULE_PATH: Search path for the module Eigen3 to be loaded by find_package
    #
    SET( EIGEN3_ROOT "./eigen-3.4.0" )
    SET( EIGEN3_DIR "${EIGEN3_ROOT}/cmake" )
    SET( EIGEN3_INCLUDE_DIR "${EIGEN3_ROOT}" )
    SET( CMAKE_MODULE_PATH EIGEN3_DIR )

    # find_package( Eigen3 3.4.0 REQUIRED )
    # include_directories is needed for the compiler to know where looking for Eigen3 header files to be included
    include_directories( ${EIGEN3_INCLUDE_DIR} )
	add_library(AliasEigen3 INTERFACE src/test_bo.cpp)
    target_include_directories(AliasEigen3 INTERFACE ${EIGEN3_INCLUDE_DIR})
    add_library(Eigen3::Eigen ALIAS AliasEigen3)

ELSE()
    find_package(Eigen3 REQUIRED)
    if((NOT TARGET Eigen3::Eigen) AND (DEFINED EIGEN3_INCLUDE_DIR))
        # Eigen 3.0, 3.1, 3.2 do not provide the target named Eigen3::Eigen
        add_library(AliasEigen3 INTERFACE src/test_bo.cpp)
        target_include_directories(AliasEigen3 INTERFACE ${EIGEN3_INCLUDE_DIR})
        add_library(Eigen3::Eigen ALIAS AliasEigen3)
    endif()
ENDIF()



# timer
set(TIMER_BUILD_TEST OFF CACHE INTERNAL "" FORCE)
add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/external/timer)

# mathtoolbox
set(MATHTOOLBOX_BUILD_EXAMPLES  OFF CACHE INTERNAL "" FORCE)
set(MATHTOOLBOX_PYTHON_BINDINGS OFF CACHE INTERNAL "" FORCE)
add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/external/mathtoolbox)

# nlopt
set(BUILD_SHARED_LIBS OFF CACHE INTERNAL "" FORCE)
set(NLOPT_GUILE       OFF CACHE INTERNAL "" FORCE)
set(NLOPT_LINK_PYTHON OFF CACHE INTERNAL "" FORCE)
set(NLOPT_MATLAB      OFF CACHE INTERNAL "" FORCE)
set(NLOPT_OCTAVE      OFF CACHE INTERNAL "" FORCE)
set(NLOPT_PYTHON      OFF CACHE INTERNAL "" FORCE)
set(NLOPT_SWIG        OFF CACHE INTERNAL "" FORCE)
set(NLOPT_TESTS       OFF CACHE INTERNAL "" FORCE)
set(NLOPT_CXX         ON  CACHE INTERNAL "" FORCE)
set(NLOPT_FORTRAN     OFF CACHE INTERNAL "" FORCE)
add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/external/nlopt)

# nlopt-util
add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/external/nlopt-util)

# parallel-util
set(PARALLELUTIL_BUILD_TESTS OFF CACHE INTERNAL "" FORCE)
add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/external/parallel-util)

###############################################
# Main library
###############################################

file(GLOB
	headers
	${CMAKE_CURRENT_SOURCE_DIR}/include/sequential-line-search/*.hpp
)
file(GLOB
	sources
	${CMAKE_CURRENT_SOURCE_DIR}/src/*.cpp
)
add_library(SequentialLineSearch STATIC ${headers} ${sources})

# for Windows
if(WIN32)
	add_definitions(-D_USE_MATH_DEFINES)
endif()

target_include_directories(SequentialLineSearch PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include)
target_link_libraries(SequentialLineSearch Eigen3::Eigen mathtoolbox nlopt nlopt-util parallel-util timer)
if(SEQUENTIAL_LINE_SEARCH_USE_NOISELESS_FORMULATION)
	target_compile_definitions(SequentialLineSearch PRIVATE SEQUENTIAL_LINE_SEARCH_USE_NOISELESS_FORMULATION)
endif()
if(SEQUENTIAL_LINE_SEARCH_USE_PARALLELIZED_MULTI_START_SEARCH)
	target_compile_definitions(SequentialLineSearch PRIVATE SEQUENTIAL_LINE_SEARCH_USE_PARALLELIZED_MULTI_START_SEARCH)
endif()

install(FILES ${headers} DESTINATION include/SequentialLineSearch)
install(TARGETS SequentialLineSearch ARCHIVE DESTINATION lib)

###############################################
# Sub-directories
###############################################

if(SEQUENTIAL_LINE_SEARCH_BUILD_COMMAND_DEMOS)
	add_subdirectory(demos/bayesian_optimization_1d)
	add_subdirectory(demos/sequential_line_search_nd)
endif()
if(SEQUENTIAL_LINE_SEARCH_BUILD_VISUAL_DEMOS)
	# Qt
	find_package(Qt5 COMPONENTS Gui Widgets Concurrent REQUIRED)

	# tinycolormap
	set(TINYCOLORMAP_WITH_EIGEN  ON  CACHE INTERNAL "" FORCE)
	set(TINYCOLORMAP_WITH_QT5    ON  CACHE INTERNAL "" FORCE)
	set(TINYCOLORMAP_BUILD_TOOLS OFF CACHE INTERNAL "" FORCE)
	add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/external/tinycolormap)

	# rand-util
	add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/external/rand-util)

	add_subdirectory(demos/bayesian_optimization_1d_gui)
	add_subdirectory(demos/bayesian_optimization_2d_gui)
	add_subdirectory(demos/preferential_bayesian_optimization_1d_gui)
	add_subdirectory(demos/sequential_line_search_2d_gui)
endif()
if(SEQUENTIAL_LINE_SEARCH_BUILD_PHOTO_DEMOS)
	# Qt
	find_package(Qt5 COMPONENTS Gui Widgets OpenGL Concurrent REQUIRED)

	# OpenGL
	find_package(OpenGL REQUIRED)

	# tinycolormap
	if(NOT TARGET tinycolormap)
		set(TINYCOLORMAP_WITH_EIGEN  ON  CACHE INTERNAL "" FORCE)
		set(TINYCOLORMAP_WITH_QT5    ON  CACHE INTERNAL "" FORCE)
		set(TINYCOLORMAP_WITH_GLM    OFF CACHE INTERNAL "" FORCE)
		set(TINYCOLORMAP_BUILD_TOOLS OFF CACHE INTERNAL "" FORCE)
		add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/external/tinycolormap)
	endif()

	# enhancer
	set(ENHANCER_USE_QT_FEATURES ON CACHE INTERNAL "" FORCE)
	set(ENHANCER_BUILD_QT_TESTS OFF CACHE INTERNAL "" FORCE)
	add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/external/enhancer)

	add_subdirectory(demos/sequential_line_search_photo)
endif()
if(SEQUENTIAL_LINE_SEARCH_BUILD_PYTHON_BINDING)
	# pybind11
	set(PYBIND11_INSTALL OFF CACHE INTERNAL "" FORCE)
	set(PYBIND11_TEST    OFF CACHE INTERNAL "" FORCE)
	add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/external/pybind11)

	add_subdirectory(python)
endif()

###############################################
# Tests
###############################################
enable_testing()
if(SEQUENTIAL_LINE_SEARCH_BUILD_COMMAND_DEMOS)
	add_test(NAME bayesian_optimization_1d_test COMMAND $<TARGET_FILE:BayesianOptimization1d>)
	add_test(NAME sequential_line_search_nd_test COMMAND $<TARGET_FILE:SequentialLineSearchNd>)
endif()


