cmake_minimum_required(VERSION 3.14)
project(QuasiplanarGraphGUI CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# CRITICAL FIX: SFML 3 requires strictly capitalized component names
find_package(SFML 3 REQUIRED COMPONENTS Graphics Window System)

add_executable(quasiplanar_gui main.cpp hds_quasiplanar.h)

# Link using the modern SFML 3 target namespaces
target_link_libraries(quasiplanar_gui PRIVATE 
    SFML::Graphics 
    SFML::Window 
    SFML::System
)
