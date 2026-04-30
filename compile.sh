#!/bin/bash
# -O3 turns on maximum performance optimizations
# -std=c++17 ensures modern C++ features are supported
# -Wall shows standard warnings to help catch bugs
g++ -O3 -Wall -std=c++17 enumeration.cpp -o ./.bin/enumeration

./.bin/enumeration
