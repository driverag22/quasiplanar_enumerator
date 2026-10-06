#!/bin/bash -l
# One task of a split search (Slurm job array); every task runs the same binary with its own task index.
# Submit from the repo root, e.g.
#   sbatch --array=0-99 --time=7-00:00:00 \
#     --export=ALL,BIN=.bin/split_c12_run1,RUN=c12_run1,TASKS=100,DEPTH=12,PASSES=1,REACH_EVERY_EDGE_FROM=1,FINAL_FACE_EVERY_EDGE_FROM=1,PLACEMENT_EVERY_EDGE_FROM=0,PLACEMENT_FILTER=1,PROGRESS_SECONDS=1800 \
#     nested/slurm/search.slurm
# Variables:
#   BIN     binary (use a frozen copy per run, all tasks and resubmissions must use the same one)
#   RUN     run name, output in $OUT/$RUN/t<task> (OUT defaults to runs)
#   TASKS   number of tasks (= array size), DEPTH split depth, PASSES 1, 2 or 12
#   PARENTS optional parents file (default: the base drawing)
#   search options (REACH_EVERY_EDGE_FROM, ...) are read from the environment by the program,
#   see nested/split_driver.h; they must be the same for all tasks of a run
#SBATCH --job-name=nested_split
#SBATCH --output=logs/%x_%A_%a.log
#SBATCH --error=logs/%x_%A_%a.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=2G
#SBATCH --partition=standard
cd "$SLURM_SUBMIT_DIR"
module load stack/2024-06 boost >/dev/null 2>&1
OUT=${OUT:-runs}
mkdir -p "$OUT/$RUN"
echo "task $SLURM_ARRAY_TASK_ID/$TASKS depth $DEPTH passes $PASSES binary $BIN host $(hostname) start $(date -Is)"
echo "options: reach_from=${REACH_EVERY_EDGE_FROM:-off} final_face_from=${FINAL_FACE_EVERY_EDGE_FROM:-off} placement_from=${PLACEMENT_EVERY_EDGE_FROM:-off} placement_filter=${PLACEMENT_FILTER:-off}"
"$BIN" search "$SLURM_ARRAY_TASK_ID" "$TASKS" "$DEPTH" "$PASSES" "$OUT/$RUN/t$SLURM_ARRAY_TASK_ID" $PARENTS
