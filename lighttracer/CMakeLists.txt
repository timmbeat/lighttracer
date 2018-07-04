cmake_minimum_required (VERSION 3.9)
project (lighttracer)

file(GLOB base_SRCS lighttracer/*.h
        lighttracer/*.cpp)


set_property(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR} PROPERTY VS_STARTUP_PROJECT lighttracer)
add_executable(lighttracer ${base_SRCS})
set_property(TARGET lighttracer PROPERTY CXX_STANDARD 17)
target_include_directories(lighttracer PRIVATE lighttracer glm/glm)