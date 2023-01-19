cmake_minimum_required ( VERSION 3.12 ) 
set ( CMAKE_CXX_STANDARD 17 ) 

project ( SurrogateModelling CXX ) 


# Add variables containing files' names 
file ( GLOB BASIS_FUNCTIONS_HEADERS library/BasisFunctions/*.hpp ) 
file ( GLOB UTILITIES_FILES library/Utilities/*.hpp ) 


# Add interface library 

add_library ( surrogatemodel INTERFACE 

    ${BASIS_FUNCTIONS_HEADERS}
    ${UTILITIES_FILES} 

)

target_include_directories ( surrogatemodel INTERFACE 

    library/BasisFunctions 
    library/Utilities 

)


# Add testrunner executable 
add_executable ( testrunner test/testrunner.cpp )
target_link_libraries ( testrunner PUBLIC surrogatemodel ) 

