cmake_minimum_required ( VERSION 3.25 ) 
set ( CMAKE_CXX_STANDARD 17 ) 
set ( CMAKE_CXX_COMPILER "g++-12" )

project ( SurrogateModelling CXX ) 

if     ( CMAKE_CXX_COMPILER_ID MATCHES "Clang" ) 
    add_compile_options ( -Wno-deprecated-declarations )
elseif ( CMAKE_CXX_COMPILER_ID MATCHES "GNU"   ) 
    add_compile_options ( -O3 -fopenmp -Wno-deprecated-declarations )
elseif ( CMAKE_CXX_COMPILER_ID MATCHES "MSVC"  )
elseif ( CMAKE_CXX_COMPILER_ID MATCHES "INTEL" )
endif ()

find_package ( Eigen3 3.4 REQUIRED NO_MODULE ) 
find_package ( OpenMP ) 
if ( OPENMP_FOUND ) 
    SET ( CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}" )
    SET ( CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}" )
endif ()

# Add interface libraries 
add_subdirectory ( library )

# Add executable 
add_executable ( main driver/main.cpp ) 
target_link_libraries ( main PUBLIC 

    basisfunctions 
    montecarlo 
    linalg  
    Eigen3::Eigen 

)

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

    enable_testing ()

    add_executable ( OverallTestrunner 

        test/testrunner.cpp 

    )

    target_link_libraries ( OverallTestrunner PUBLIC 

        basisfunctions 
        montecarlo 
        linalg  
        GTest::gtest_main 

    ) 

    include ( GoogleTest ) 
    gtest_discover_tests ( OverallTestrunner ) 

endif ()

