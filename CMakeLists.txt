cmake_minimum_required(VERSION 2.6)

SET(ENV{BOOST_ROOT} "/usr/include/boost")

ADD_DEFINITIONS(
    -std=c++11
    # Other flags
)
find_package(Boost EXACT REQUIRED COMPONENTS python3)

INCLUDE_DIRECTORIES("${Boost_INCLUDE_DIRS}" "/usr/include/python3.4m/" "./")

ADD_LIBRARY(mathieu SHARED mathieu.cpp mathieu.h mathieu_wrapper.h mathieu_wrapper.cpp main.cpp)
TARGET_LINK_LIBRARIES(mathieu ${Boost_LIBRARIES})
