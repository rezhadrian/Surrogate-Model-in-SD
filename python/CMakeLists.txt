cmake_minimum_required ( VERSION 3.25 ) 
project ( SMSD_python CXX ) 


# ----- Set up library for C++ bindings 
set( PYTHON_BINDINGS_FOLDER bindings )


# ----- Find all C++ files in bindings folder
file( GLOB BINDING_SOURCES ${PYTHON_BINDINGS_FOLDER}/*.cpp )
file( GLOB BINDING_HEADERS ${PYTHON_BINDINGS_FOLDER}/*.h* )


# ----- Add shared library of python bindings
pybind11_add_module( SMSD ${BINDING_SOURCES} ${BINDING_HEADERS} )

target_link_libraries( SMSD PRIVATE surrogatemodel )


# ----- Specify relative path the shared library shall be installed to
install( TARGETS SMSD LIBRARY DESTINATION ${CMAKE_INSTALL_PREFIX} )
  

# ----- Set up copying python scipts 
set( PYTHON_SCRIPT_FOLDER scripts )

# ----- Find all python files in script folder
file( GLOB PYTHON_SCRIPTS ${PYTHON_SCRIPT_FOLDER}/*.py )

# ----- Add all files to install target
foreach( filename ${PYTHON_SCRIPTS} )
  install(FILES ${filename} DESTINATION ${CMAKE_INSTALL_PREFIX} )
endforeach()


