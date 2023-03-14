cmake_minimum_required ( VERSION 3.25 ) 
project ( SurrogateModel CXX ) 

set ( CMAKE_CXX_STANDARD 17 ) 


# ----- Add package "Eigen3" and "Boost" 

find_package ( Eigen3 3.4 REQUIRED NO_MODULE ) 
find_package ( Boost 1.81 REQUIRED ) 

# find_package ( OpenMP ) 
# if ( OPENMP_FOUND ) 
#     SET ( CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}" )
#     SET ( CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}" )
# endif ()


# ----- Add "modules" subdirectory, contains templated functionalities 

add_subdirectory ( modules )


# ----- Add shared library "surrogatemodel" 

add_library ( surrogatemodel SHARED 

    declarations/Surrogate_MassSpringDamper.hpp 

    implementations/MassSpringDamper_DMCS.cpp 
    implementations/MassSpringDamper_IPCE.cpp 
    implementations/MassSpringDamper_IRPCE.cpp 
    implementations/MassSpringDamper_NIPCE.cpp 
    implementations/MassSpringDamper_NIRPCE.cpp 

)

target_include_directories ( surrogatemodel PUBLIC 

    ${CMAKE_CURRENT_SOURCE_DIR}/declarations 
    ${CMAKE_CURRENT_SOURCE_DIR}/utility  
    ${Boost_INCLUDE_DIR} 

)

target_link_libraries ( surrogatemodel PUBLIC 

    surrogatemodule 
    analyticalmodel 
    Eigen3::Eigen 
    ${Boost_LIBRARIES} 

)

# specify the relative path the shared library object shall be installed to
if( WIN32 )
  install( TARGETS surrogatemodel RUNTIME DESTINATION ${CMAKE_INSTALL_PREFIX} )
else( )
  install( TARGETS surrogatemodel LIBRARY DESTINATION ${CMAKE_INSTALL_PREFIX} )
endif( )


# ----- Add executable 

# add_executable ( main driver/main.cpp ) 
# target_link_libraries ( main PUBLIC 
#
#     surrogatemodel
#     Eigen3::Eigen 
#
# )

# Unit Test using Google Test 
# option ( PROJECT_TEST "Enable Google Test for whole project" OFF )
#
# if ( PROJECT_TEST ) 
#
#     if ( NOT SMSD_GTEST )
#
#         include ( FetchContent )
#
#         FetchContent_Declare(
#             googletest
#             URL https://github.com/google/googletest/archive/03597a01ee50ed33e9dfd640b249b4be3799d395.zip
#         )
#         
#         set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)
#         FetchContent_MakeAvailable ( googletest ) 
#        
#     endif ()
#
#     set ( SMSD_GTEST ON )
#
#     enable_testing ()
#
#     add_executable ( OverallTestrunner 
#
#         test/testrunner.cpp 
#
#     )
#
#     target_link_libraries ( OverallTestrunner PUBLIC 
#
#         basisfunctions 
#         montecarlo 
#         GTest::gtest_main 
#
#     ) 
#
#     include ( GoogleTest ) 
#     gtest_discover_tests ( OverallTestrunner ) 
#
# endif ()

