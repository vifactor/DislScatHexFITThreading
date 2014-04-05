cmake_minimum_required(VERSION 2.8)

project(FitSkew)

#EXTERNAL LIBRARIES
include(ExternalProject)

#port library
ExternalProject_Add(
    port_project
    PREFIX ${CMAKE_BINARY_DIR}
    URL ${CMAKE_CURRENT_SOURCE_DIR}/tarballs/port.tar.gz
)

#levmar library
ExternalProject_Add(
    levmar_project
    PREFIX ${CMAKE_BINARY_DIR}
    URL http://users.ics.forth.gr/~lourakis/levmar/levmar-2.6.tgz
    INSTALL_COMMAND ""
)

#config library
set(config_project_prefix "${CMAKE_BINARY_DIR}")
ExternalProject_Add(
    config_project
    PREFIX ${config_project_prefix}
    URL http://www.hyperrealm.com/libconfig/libconfig-1.4.9.tar.gz
    SOURCE_DIR ${config_project_prefix}/src/config_project
    CONFIGURE_COMMAND ${config_project_prefix}/src/config_project/configure --prefix=${config_project_prefix}/src/config_project-build
    BUILD_COMMAND make
    BUILD_IN_SOURCE 1
)

#compile project itself
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/src)
include_directories(${CMAKE_BINARY_DIR}/src/port_project-build)
include_directories(${CMAKE_BINARY_DIR}/src/config_project-build/include)
include_directories(${CMAKE_BINARY_DIR}/src/levmar_project)

link_directories(${CMAKE_BINARY_DIR}/src/port_project-build)
link_directories(${CMAKE_BINARY_DIR}/src/config_project-build/lib)
link_directories(${CMAKE_BINARY_DIR}/src/levmar_project-build)

file(GLOB CPP_SOURCES "${CMAKE_CURRENT_SOURCE_DIR}/src/*.cpp")
add_executable(FitSkew ${CPP_SOURCES})

add_dependencies(FitSkew port_project)
add_dependencies(FitSkew levmar_project)
add_dependencies(FitSkew config_project)

target_link_libraries(FitSkew port levmar config++ lapack blas gsl)
