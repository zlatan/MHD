cmake_minimum_required(VERSION 3.10)
project(mhd_disk C)

set(CMAKE_C_STANDARD 11)

add_executable(mhd_disk main.c)
target_link_libraries(mhd_disk m)