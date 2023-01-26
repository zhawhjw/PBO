cmake_minimum_required(VERSION 3.1)
install(FILES ${CMAKE_CURRENT_SOURCE_DIR}/include/nlopt-util.hpp DESTINATION include)

IF (WIN32)
    SET( EIGEN3_ROOT "../../eigen-3.4.0" )
    SET( EIGEN3_DIR "${EIGEN3_ROOT}/cmake" )
    SET( EIGEN3_INCLUDE_DIR "${EIGEN3_ROOT}" )
    SET( CMAKE_MODULE_PATH EIGEN3_DIR )
    include_directories( ${EIGEN3_INCLUDE_DIR} )

ELSE()
    find_package(Eigen3)
    
ENDIF()

add_library(nlopt-util INTERFACE)
target_sources(nlopt-util INTERFACE ${CMAKE_CURRENT_SOURCE_DIR}/include/nlopt-util.hpp)
target_include_directories(nlopt-util INTERFACE ${CMAKE_CURRENT_SOURCE_DIR}/include)
target_link_libraries(nlopt-util INTERFACE Eigen3::Eigen)




