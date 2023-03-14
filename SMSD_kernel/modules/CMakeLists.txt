cmake_minimum_required ( VERSION 3.25 ) 
set ( CMAKE_CXX_STANDARD 17 ) 

project ( SurrogateModule CXX ) 

# Add interface libraries 
add_subdirectory ( AnalyticalModel ) 
add_subdirectory ( BasisFunctions )
add_subdirectory ( MonteCarlo )

# Combine interface libraries 
add_library ( surrogatemodule INTERFACE )
target_link_libraries ( surrogatemodule INTERFACE 

    analyticalmodel 
    basisfunctions 
    montecarlo

)

# Unit Test using Google Test 
option ( SURROGATE_TEST "Enable Google Test for all libraries" ON )

if ( SURROGATE_TEST ) 

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

endif ()
