cmake_minimum_required ( VERSION 3.25 ) 
set ( CMAKE_CXX_STANDARD 17 ) 

project ( BasisFunctions CXX ) 

# Add interface library "basisfunctions" 
add_library ( basisfunctions INTERFACE 

    declarations/BasisFunctions.hpp 

    implementations/HermitePolynomials_imp.hpp 
    implementations/MultiIndex_imp.hpp 
    implementations/SupplementaryMaths_imp_BF.hpp 
    implementations/Truncations_imp.hpp 

    utility/LibrariesLoader_BF.hpp

)

target_include_directories ( basisfunctions INTERFACE 

    ${CMAKE_CURRENT_SOURCE_DIR}/declarations 
    ${CMAKE_CURRENT_SOURCE_DIR}/implementations  
    ${CMAKE_CURRENT_SOURCE_DIR}/utility

)


# Unit Test using Google Test 
option ( BASIS_FUNCTION_TEST "Enable Google Test for Basis Functions" ON )

if ( BASIS_FUNCTION_TEST ) 

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

    add_executable ( BasisFunctions_testrunner 

        test/HermitePolynomial_test.cpp 
        test/HermitePolynomials_test.cpp 
        test/MultiIndex_test.cpp 
        test/SupplementaryMaths_BF_test.cpp 
        test/TotalTruncations_test.cpp 

    )

    target_include_directories ( BasisFunctions_testrunner PUBLIC 

        ${CMAKE_CURRENT_SOURCE_DIR}/test 

    )

    target_link_libraries ( BasisFunctions_testrunner PUBLIC 

        basisfunctions 
        GTest::gtest_main 

    ) 

    include ( GoogleTest ) 
    gtest_discover_tests ( BasisFunctions_testrunner ) 
 
endif ()
