cmake_minimum_required ( VERSION 3.25 ) 
set ( CMAKE_CXX_STANDARD 17 ) 

project ( MonteCarlo CXX ) 

find_package ( Eigen3 3.4 REQUIRED NO_MODULE ) 
find_package ( Boost 1.81 REQUIRED )

if     ( CMAKE_CXX_COMPILER_ID MATCHES "Clang" ) 
    add_compile_options ( -Wno-deprecated-declarations )
elseif ( CMAKE_CXX_COMPILER_ID MATCHES "GNU"   ) 
    add_compile_options ( -Wno-deprecated-declarations )
elseif ( CMAKE_CXX_COMPILER_ID MATCHES "MSVC"  )
elseif ( CMAKE_CXX_COMPILER_ID MATCHES "INTEL" )
endif ()

# Add interface library "montecarlo" 
add_library ( montecarlo INTERFACE 

    declarations/MonteCarlo.hpp 

    declarations/include/StatisticDist.hpp 
    declarations/include/StatisticDist_imp.hpp 

    implementations/AnalyticalModel_imp.hpp 
    implementations/LatinHypercubeSampling_imp.hpp 
    implementations/VariableGeneration_imp.hpp 

)

target_include_directories ( montecarlo INTERFACE 

    ${CMAKE_CURRENT_SOURCE_DIR}/declarations 
    ${CMAKE_CURRENT_SOURCE_DIR}/declarations/include 
    ${CMAKE_CURRENT_SOURCE_DIR}/implementations 
    ${CMAKE_CURRENT_SOURCE_DIR}/utility   

    ${Boost_INCLUDE_DIR} 

)

target_link_libraries ( montecarlo INTERFACE 

    Eigen3::Eigen 
    ${Boost_LIBRARIES} 

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

        test/AnalyticalModel_test.cpp 
        test/LatinHypercubeSampling_test.cpp
        test/VariableGeneration_test.cpp 

        test/include/StatisticDist_test.cpp 

    )

    target_link_libraries ( MonteCarlo_testrunner PUBLIC 

        montecarlo
        GTest::gtest_main 

    ) 

    include ( GoogleTest ) 
    gtest_discover_tests ( MonteCarlo_testrunner ) 
 
endif ()
