cmake_minimum_required ( VERSION 3.25 ) 
set ( CMAKE_CXX_STANDARD 17 ) 

project ( MonteCarlo CXX ) 

# Add interface library "montecarlo" 
add_library ( montecarlo INTERFACE 

    declarations/InverseCDF.hpp 
    declarations/LatinHypercubeSampling.hpp 
    declarations/SpecialMatrix.hpp 
    declarations/SupplementaryMaths_MC.hpp
    implementations/InverseCDF_imp.hpp 
    implementations/LatinHypercubeSampling_imp.hpp 
    implementations/SpecialMatrix_imp.hpp 
    implementations/SupplementaryMaths_imp_MC.hpp 

)

target_include_directories ( montecarlo INTERFACE 

    ${CMAKE_CURRENT_SOURCE_DIR}/declarations 
    ${CMAKE_CURRENT_SOURCE_DIR}/implementations 
    ${CMAKE_CURRENT_SOURCE_DIR}/utility   

)


# Unit Test using Google Test 
option ( MONTE_CARLO_TEST "Enable Google Test for Monte Carlo" ON )

if ( MONTE_CARLO_TEST ) 

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

    add_executable ( MonteCarlo_testrunner

        test/CholeskyDecomposition_test.cpp 
        test/InverseCDF_test.cpp 
        test/LatinHypercubeSampling_test.cpp
        test/SpecialMatrix_test.cpp 

    )

    target_include_directories ( MonteCarlo_testrunner PUBLIC 

        ${CMAKE_CURRENT_SOURCE_DIR}/test 

    )

    target_link_libraries ( MonteCarlo_testrunner PUBLIC 

        montecarlo
        GTest::gtest_main 

    ) 

    include ( GoogleTest ) 
    gtest_discover_tests ( MonteCarlo_testrunner ) 
 
endif ()
