#!/bin/bash
# Build all validation drivers (run from anywhere). Binaries go to validation/.bin/ (at most 6 compiles at once)
set -e
cd "$(dirname "$0")/.."
F="-O3 -Wall -std=c++17 -Iinclude"
mkdir -p validation/.bin
{
  for k in 1 2 3 4 5; do echo "-DKLIM=$k validation/enum.cpp -o validation/.bin/enum_k$k"; done
  # quasiplane: klim = 2n-7 for n = 4,...,9 vertices (enum_q<klim>)
  for k in 1 3 5 7 9 11; do echo "-DQUASI -DKLIM=$k validation/enum.cpp -o validation/.bin/enum_q$k"; done
  [ -f validation/hrw.cpp ] && echo "validation/hrw.cpp -o validation/.bin/hrw"
  [ -f validation/pruning.cpp ] && echo "validation/pruning.cpp -o validation/.bin/pruning"
} | xargs -P 6 -I{} sh -c "g++ $F {}"
echo built
