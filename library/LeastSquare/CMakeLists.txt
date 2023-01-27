cmake_minimum_required ( VERSION 3.25 ) 
set ( CMAKE_CXX_STANDARD 17 ) 

project ( LeastSquare CXX ) 

add_subdirectory ( dependencies ) 

# Add interface library "montecarlo" 
add_library ( leastsquare INTERFACE 

    declarations/LeastSquare.hpp 

)

target_include_directories ( leastsquare INTERFACE 

    ${CMAKE_CURRENT_SOURCE_DIR}/declarations 
    ${CMAKE_CURRENT_SOURCE_DIR}/declarations/include 
    ${CMAKE_CURRENT_SOURCE_DIR}/implementations 
    ${CMAKE_CURRENT_SOURCE_DIR}/utility   

)

target_link_libraries ( leastsquare INTERFACE 

    linalg

)


# Unit Test using Google Test 
option ( LEAST_SQUARE_TEST "Enable Google Test for Least Square" ON )

if ( LEAST_SQUARE_TEST ) 

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

    add_executable ( LeastSquare_testrunner

        test/LeastSquare_test.cpp 

    )

    target_link_libraries ( LeastSquare_testrunner PUBLIC 

        leastsquare
        GTest::gtest_main 

    ) 

    include ( GoogleTest ) 
    gtest_discover_tests ( LeastSquare_testrunner ) 
 
endif ()
