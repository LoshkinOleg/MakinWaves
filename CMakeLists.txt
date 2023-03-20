cmake_minimum_required (VERSION 3.14)
project(WaveSim2d)

# Using C++20 .
set (CMAKE_CXX_STANDARD 20)
set (CMAKE_CXX_STANDARD_REQUIRED ON)

# Make an out-of-source CMake build under /build for relative paths to work as intended.
if(NOT EXISTS "${PROJECT_SOURCE_DIR}/build/")
	message(FATAL_ERROR "Please specify an out-of-source directory ‘build/’ in the project’s root directory. If you don’t know what an out-of-source build is, here’s a link: https://cmake.org/cmake/help/book/mastering-cmake/chapter/Getting%20Started.html")
endif()

# Detect platform.
if (APPLE)
	message(FATAL_ERROR "Scout's development on Apple platforms isn't supported.")
endif()
if (UNIX AND NOT APPLE)
	message(FATAL_ERROR "Scout's development on Linux platforms isn't supported.")
endif()


# Thirdparty headers locations.
set(RayLib_IncludePath ${PROJECT_SOURCE_DIR}/thirdparty/raylib/include/)

# Third party prebuild static libraries.
if(WIN32)

set(RayLib_Libs
	${PROJECT_SOURCE_DIR}/thirdparty/raylib/lib/windows/raylib.lib
	)
	
endif()

# Application's headers location.
set(Application_IncludePath ${PROJECT_SOURCE_DIR}/application/include/)

# Application's sources.
file(GLOB Application_Inc ${PROJECT_SOURCE_DIR}/application/include/*.h)
file(GLOB Application_Src ${PROJECT_SOURCE_DIR}/application/src/*.cpp)

# Application's executable declaration.
add_executable(Application ${Application_Inc} ${Application_Src})
set_target_properties(Application PROPERTIES LINKER_LANGUAGE CXX)

# Application's include directories.
target_include_directories(Application PUBLIC
	${Application_IncludePath}
	${RayLib_IncludePath}
)

# Linking specifications for Application.
if(WIN32)

target_link_libraries(Application PUBLIC
	general winmm.lib # A dependancy of raylib for windows.
	general ${RayLib_Libs}
)

endif()