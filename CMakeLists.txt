cmake_minimum_required (VERSION 2.6)
project (RobustDubins)

# set flags
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS}  -std=c++11 -fPIC")
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/lib)

# get path variables
get_filename_component(ROOT_DIR ${CMAKE_CURRENT_SOURCE_DIR} DIRECTORY)

# include headers into the project
include_directories(
${CMAKE_CURRENT_SOURCE_DIR}/src
# custom packages
# ${ROOT_DIR}/MathTools/src
)

# external libraries
link_directories(
#/lib
#/usr/lib
#/usr/local/lib
# custom packages
# ${ROOT_DIR}/MathTools/lib
${ROOT_DIR}/RobustDubins/lib
)


# use function file(GLOB...) to grab all source files using wildcard syntax
file(GLOB SOURCES "${CMAKE_CURRENT_SOURCE_DIR}/src/*.cpp")
add_library(RobustDubins ${SOURCES})
target_link_libraries(RobustDubins) 



# build programs
  file(GLOB PROGRAMS "${CMAKE_CURRENT_SOURCE_DIR}/programs/*.cpp")  
  foreach( testsourcefile ${PROGRAMS} )
      # I used a simple string replace, to cut off .cpp.
      get_filename_component(testname ${testsourcefile} NAME_WE)
      add_executable(${testname} ${testsourcefile} )
      set_target_properties(${testname} PROPERTIES RUNTIME_OUTPUT_DIRECTORY "../bin")
      # Make sure YourLib is linked to each app
      target_link_libraries(${testname} RobustDubins)
  endforeach()

