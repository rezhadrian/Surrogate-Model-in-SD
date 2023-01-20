cmake_minimum_required ( VERSION 3.25 ) 
set ( CMAKE_CXX_STANDARD 17 ) 

project ( SurrogateModelling CXX ) 

# Add interface libraries 
add_subdirectory ( library )

# Unit Test using Google Test 
option ( PROJECT_TEST "Enable Google Test for whole project" ON )

if ( PROJECT_TEST ) 

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

