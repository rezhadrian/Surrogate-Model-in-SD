cmake_minimum_required ( VERSION 3.25 ) 
set ( CMAKE_CXX_STANDARD 17 ) 

project ( AnalyticalModel CXX ) 

find_package ( Eigen3 3.4 REQUIRED NO_MODULE ) 
find_package ( Boost 1.81 REQUIRED )

if     ( CMAKE_CXX_COMPILER_ID MATCHES "Clang" ) 
    add_compile_options ( -Wno-deprecated-declarations )
elseif ( CMAKE_CXX_COMPILER_ID MATCHES "GNU"   ) 
    add_compile_options ( -Wno-deprecated-declarations )
elseif ( CMAKE_CXX_COMPILER_ID MATCHES "MSVC"  )
elseif ( CMAKE_CXX_COMPILER_ID MATCHES "INTEL" )
endif ()

# Add interface library "basisfunctions" 
add_library ( analyticalmodel INTERFACE 

    declarations/AnalyticalModel.hpp 

    implementations/MassSpringDamper_imp.hpp 

    utility/LibrariesLoader_AM.hpp

)

target_include_directories ( analyticalmodel INTERFACE 

    ${CMAKE_CURRENT_SOURCE_DIR}/declarations 
    ${CMAKE_CURRENT_SOURCE_DIR}/implementations  
    ${CMAKE_CURRENT_SOURCE_DIR}/utility

    ${Boost_INCLUDE_DIR} 

)

target_link_libraries ( analyticalmodel INTERFACE 

    Eigen3::Eigen 
    ${Boost_LIBRARIES} 

)


# Unit Test using Google Test 
option ( ANALYTICAL_MODEL_TEST "Enable Google Test for Basis Functions" ON )

if ( ANALYTICAL_MODEL_TEST ) 

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

    add_executable ( AnalyticalModel_testrunner 

        test/MassSpringDamper_test.cpp 

    )

    target_link_libraries ( AnalyticalModel_testrunner PUBLIC 

        analyticalmodel
        GTest::gtest_main 

    ) 

    include ( GoogleTest ) 
    gtest_discover_tests ( AnalyticalModel_testrunner ) 
 
endif ()
