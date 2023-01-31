cmake_minimum_required ( VERSION 3.25 ) 
set ( CMAKE_CXX_STANDARD 17 ) 

project ( LinearAlgebra CXX ) 

# Add interface library "linalg" 
add_library ( linalg INTERFACE 

    declarations/LinearAlgebra.hpp 

    declarations/classes/Matrix.hpp 
    declarations/classes/operations/Arithmetic.hpp 
    declarations/classes/operations/InnerProducts.hpp 
    declarations/classes/operations/Norms.hpp 

    declarations/decompositions/Cholesky.hpp 

    declarations/solvers/GaussSeidel.hpp 

    implementations/classes/Matrix_imp.hpp 
    implementations/classes/operations/Arithmetic_imp.hpp 
    implementations/classes/operations/InnerProducts_imp.hpp 
    implementations/classes/operations/Norms_imp.hpp 

    implementations/decompositions/Cholesky_imp.hpp 

    implementations/solvers/GaussSeidel_imp.hpp 


)

target_include_directories ( linalg INTERFACE 

    ${CMAKE_CURRENT_SOURCE_DIR}/declarations 
    ${CMAKE_CURRENT_SOURCE_DIR}/declarations/classes 
    ${CMAKE_CURRENT_SOURCE_DIR}/declarations/classes/operations 

    ${CMAKE_CURRENT_SOURCE_DIR}/declarations/decompositions 

    ${CMAKE_CURRENT_SOURCE_DIR}/declarations/solvers 
 

    ${CMAKE_CURRENT_SOURCE_DIR}/implementations 
    ${CMAKE_CURRENT_SOURCE_DIR}/implementations/classes 
    ${CMAKE_CURRENT_SOURCE_DIR}/implementations/classes/operations

    ${CMAKE_CURRENT_SOURCE_DIR}/implementations/decompositions 

    ${CMAKE_CURRENT_SOURCE_DIR}/implementations/solvers 

    ${CMAKE_CURRENT_SOURCE_DIR}/utility   

)


# Unit Test using Google Test 
option ( LINEAR_ALGEBRA_TEST "Enable Google Test for Least Square" ON )

if ( LINEAR_ALGEBRA_TEST ) 

    if ( NOT SMSD_GTEST )

        include ( FetchContent )

        FetchContent_Declare(
            googletest
            URL https://github.com/google/googletest/archive/03597a01ee50ed33e9dfd640b249b4be3799d395.zip
        )

        set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)
        FetchContent_MakeAvailable ( googletest ) 
       
    endif ()

    set ( SMSD_GTEST ON )

    enable_testing ()

    add_executable ( LinAlg_testrunner

        test/Cholesky_test.cpp 
        test/GaussSeidel_test.cpp 
        test/LeastSquare_test.cpp 

    )

    target_link_libraries ( LinAlg_testrunner PUBLIC 

        linalg
        GTest::gtest_main 

    ) 

    include ( GoogleTest ) 
    gtest_discover_tests ( LinAlg_testrunner ) 
 
endif ()
