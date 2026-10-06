#!/bin/bash -l
# Build the split drivers (submit from the repo root: sbatch nested/slurm/build.slurm).
# Compiling on the login node works as well (bash -l; module load stack/2024-06 boost; same g++ line).
# CYCLES: which drivers to build (default "12 14"); output: .bin/split_c<k>
#SBATCH --job-name=nested_build
#SBATCH --output=logs/%x_%j.log
#SBATCH --error=logs/%x_%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=2G
#SBATCH --time=00:20:00
#SBATCH --partition=standard
cd "$SLURM_SUBMIT_DIR"
module load stack/2024-06 boost >/dev/null 2>&1
mkdir -p .bin
for C in ${CYCLES:-12 14}; do
    g++ -O3 -Wall -std=c++17 -Iinclude nested/iterative_nested_c${C}_split.cpp -o .bin/split_c${C}
done
