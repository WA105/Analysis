project(Analysis)
cmake_minimum_required(VERSION 2.8)

#Set the C++ standard to c++11.
#set(CMAKE_CXX_STANDARD 11)

#Enable the use of runtime paths on Mac OSX.
#set(CMAKE_MACOSX_RPATH ON)

#ROOT libraries
list(APPEND CMAKE_PREFIX_PATH $ENV{ROOTSYS})
find_package(ROOT REQUIRED)
include(${ROOT_USE_FILE})
include_directories(${ROOT_INCLUDE_DIRS})
link_directories(${ROOT_LIBRARY_DIR})

#Bring the headers
include_directories(../include)

#Add the CMakeLists.txt on the subdirectory
add_subdirectory(../source ${PROJECT_SOURCE_DIR}/build/sourcebuild)

#Compile custom libraries
file(GLOB sources ${PROJECT_SOURCE_DIR}/../source/*.cc)
file(GLOB headers ${PROJECT_SOURCE_DIR}/../include/*.hh)

#add here the name of your executable
add_executable(noise ${PROJECT_SOURCE_DIR}/noise.cc ${sources} ${headers})
target_link_libraries(noise ${ROOT_LIBRARIES} sqlite3 CustomClasses)
