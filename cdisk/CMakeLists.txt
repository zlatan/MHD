cmake_minimum_required(VERSION 3.10)
project(cdisk C)

set(CMAKE_C_STANDARD 11)

add_executable(cdisk main.c)
target_link_libraries(cdisk m)
