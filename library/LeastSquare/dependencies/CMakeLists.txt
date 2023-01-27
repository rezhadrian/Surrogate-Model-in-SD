cmake_minimum_required ( VERSION 3.25 ) 
set ( CMAKE_CXX_STANDARD 17 ) 

project ( LinAlg CXX ) 

# Add interface library "montecarlo" 
add_library ( linalg INTERFACE 

    classes/Matrix.hpp 
    classes/Matrix_imp.hpp 
    classes/operations/Arithmetic.hpp 
    classes/operations/Arithmetic_imp.hpp 
    classes/operations/Norms.hpp 
    classes/operations/Norms_imp.hpp 
    solvers/GaussSeidel.hpp 
    solvers/GaussSeidel_imp.hpp 
    utility/LibrariesLoader_LA.hpp 

)

target_include_directories ( linalg INTERFACE 

    ${CMAKE_CURRENT_SOURCE_DIR}/classes   
    ${CMAKE_CURRENT_SOURCE_DIR}/classes/operations 
    ${CMAKE_CURRENT_SOURCE_DIR}/solvers 
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

        test/GaussSeidel_test.cpp 

    )

    target_link_libraries ( LinAlg_testrunner PUBLIC 

        linalg
        GTest::gtest_main 

    ) 

    include ( GoogleTest ) 
    gtest_discover_tests ( LinAlg_testrunner ) 
 
endif ()
