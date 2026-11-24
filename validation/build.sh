#!/bin/bash
# Builds the validation program enum.cpp for every crossing bound used in the runs (run from anywhere).
#   enum_k<k>: simple k-plane drawings, k = 1..5
#   enum_q<k>: simple quasiplane drawings with at most k crossings per edge; k = 2n-7 for K_n (Lemma 1)
# Binaries go to validation/.bin/. The runs are listed in logs/q_*.txt, one command per line, to be run from
# validation/, e.g.  while read -r c; do echo "=== $c"; $c; done < logs/q_small.txt > logs/small.log
set -e
cd "$(dirname "$0")/.."
F="-O3 -std=c++17 -Iinclude -pthread"
mkdir -p validation/.bin
{
  for k in 1 2 3 4 5; do echo "-DKLIM=$k validation/enum.cpp -o validation/.bin/enum_k$k"; done
  for k in 1 3 5 7 9 11; do echo "-DQUASI -DKLIM=$k validation/enum.cpp -o validation/.bin/enum_q$k"; done
} | xargs -P 6 -I{} sh -c "g++ $F {}"
echo built
