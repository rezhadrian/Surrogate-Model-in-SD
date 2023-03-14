cmake_minimum_required ( VERSION 3.25 ) 
project ( SMSD CXX ) 

set ( CMAKE_CXX_STANDARD 17 ) 


# ----- Add compile options based on compiler id 

if ( CMAKE_CXX_COMPILER_ID MATCHES "Clang" ) 
    add_compile_options ( -Wno-deprecated-declarations ) 
elseif ( CMAKE_CXX_COMPILER_ID MATCHES "GNU" ) 
    add_compile_options ( -Wno-deprecated-declarations ) 
elseif ( CMAKE_CXX_COMPILER_ID MATCHES "MSVC" ) 
elseif ( CMAKE_CXX_COMPILER_ID MATCHES "INTEL" ) 
endif ()


# ----- Set up install path 

if ( CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT )

    set ( 
            CMAKE_INSTALL_PREFIX 
            "${CMAKE_BINARY_DIR}/install" 
            CACHE 
            PATH 
            "default install path" 
            FORCE 
    )

endif( )

SET( CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}" )


# ----- Add "SMSD_kernel" subdirectory, contains main library 

add_subdirectory ( SMSD_kernel ) 


# ----- Set python bindings 

set(PYBIND11_CPP_STANDARD -std=c++17)
set(PYBIND11_INSTALL ON CACHE BOOL "Enable pybind11 install" )

add_subdirectory( external/pybind11 )
add_subdirectory( python )

