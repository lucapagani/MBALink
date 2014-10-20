cmake_minimum_required(VERSION 2.6)
project(mba_link)

set ( CMAKE_BUILD_TYPE "Debug" )

add_definitions("-std=c++11 -Wall -Wno-sign-compare")

include_directories(/usr/local/include/mba /usr/local/include/lsmg/include)

include_directories(./include)

add_library ( mba_link ./src/mba_link.cpp ./src/util_mba_link.cpp )

# add_executable(mba_link)

# install(TARGETS mba_link RUNTIME DESTINATION bin)

install(TARGETS mba_link
        DESTINATION "/usr/local/lib"
        )
INSTALL(FILES ./include/mba_link.hpp ./include/util_mba_link.hpp
        DESTINATION "/usr/local/include"
        )
