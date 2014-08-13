cmake_minimum_required(VERSION 2.8.4)
project(spherepack C)

if (NOT ${CMAKE_CURRENT_SOURCE_DIR} STREQUAL ${CMAKE_SOURCE_DIR})
    set(SPHEREPACK_INCLUDE ${CMAKE_CURRENT_SOURCE_DIR}/src PARENT_SCOPE)
endif ()

include_directories(contrib)
add_subdirectory(contrib)
add_subdirectory(src)

