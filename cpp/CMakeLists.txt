cmake_minimum_required(VERSION 3.4)
project(graph_algorithms)

set(CMAKE_CXX_STANDARD 17)

include_directories(/home/fanta/.local/boost_1_70_0 /home/fanta/.local/include)
link_directories(/home/fanta/.local/boost_1_70_0/stage/lib)

if (MSVC)
    add_definitions(-D_CRT_SECURE_NO_WARNINGS -wd4018 -wd4244 -wd4800 -wd4805)
else()
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Werror -Wno-unused-variable -Wno-sign-compare")
    link_libraries(pthread)
endif()

set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -O0")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O0")

add_executable(graphs graphs.cpp)

set(CMAKE_CXX_LINK_EXECUTABLE "${CMAKE_CXX_LINK_EXECUTABLE} -lstdc++fs -lcurl -lboost_iostreams")
