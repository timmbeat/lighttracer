cmake_minimum_required (VERSION 2.8)
project (lighttracer)

add_subdirectory(glm)
#Headers
include_directories(lighttracer)

set(SOURCES lighttracer/mcml_parser.cpp lighttracer/propagation.cpp lighttracer/stdafx.cpp)

add_executable(lighttracer/propagation.cpp)

file(GLOB base_SRCS lighttracer/*.h
		    lighttracer/*.cpp)

ADD_EXECUTABLE(${PROJECT_NAME} ${base_SRCS})

