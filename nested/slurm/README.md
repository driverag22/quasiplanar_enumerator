# Running the split nested-cycle search on a Slurm cluster (Euler)

All commands are run from the **repo root**; the scripts `cd` to the directory `sbatch` was called from.
Slurm writes its logs to `logs/` there, which has to exist: `mkdir -p logs`.
Code: `nested/iterative_nested_c12_split.cpp` / `_c14_split.cpp` (gadget), `nested/iterative_nested_split.h` (search),
`nested/split_driver.h` (command line: `search`, `count`, `merge`, `merge_solutions`).

## 1. Build
On the login node (or `sbatch nested/slurm/build.slurm`):

    bash -l
    module load stack/2024-06 boost
    mkdir -p .bin
    g++ -O3 -Wall -std=c++17 -Iinclude nested/iterative_nested_c12_split.cpp -o .bin/split_c12

Freeze a copy per run (`cp .bin/split_c12 .bin/split_c12_run1`) and use it for all tasks, resubmissions and the
post job of that run: the prefix numbering must be identical in every task.

## 2. Choose the split depth
The options are read from the environment; set them exactly as for the run:

    export REACH_EVERY_EDGE_FROM=1 FINAL_FACE_EVERY_EDGE_FROM=1 PLACEMENT_EVERY_EDGE_FROM=0 PLACEMENT_FILTER=1
    .bin/split_c12_run1 count 12 1        # prints the number of prefixes at depth 12 (= EXPECTED below)

With these options: C12 depth 12 -> 876,629 prefixes, C14 depth 12 -> 877,520.

## 3. Pass 1 on the base drawing (job array)

    sbatch --array=0-99 --time=7-00:00:00 \
      --export=ALL,BIN=.bin/split_c12_run1,RUN=c12_run1,TASKS=100,DEPTH=12,PASSES=1,REACH_EVERY_EDGE_FROM=1,FINAL_FACE_EVERY_EDGE_FROM=1,PLACEMENT_EVERY_EDGE_FROM=0,PLACEMENT_FILTER=1,PROGRESS_SECONDS=1800 \
      nested/slurm/search.slurm

Task i searches below the prefixes j with j mod TASKS == i and writes `runs/c12_run1/t<i>/`:
`children.jsonl` (intermediate drawings, deduplicated per task only), `done_units.txt` (finished prefixes),
`extensible.txt`. If tasks are killed (time limit), resubmit the same command: finished prefixes are skipped.

## 4. Merge + Pass 2 (one job, after the array)

    sbatch --dependency=afterany:<array job id> \
      --export=ALL,BIN=.bin/split_c12_run1,RUN=c12_run1,EXPECTED=876629 nested/slurm/post.slurm

Checks that all prefixes are done, then writes
- `runs/c12_run1_merged/parents.jsonl`: distinct intermediate drawings (gadget-aware isomorphism),
- `runs/c12_run1_final_2cycles/`: Pass 2 on the base drawing, deduplicated (strong isomorphism),
- `runs/c12_run1_final_3cycles/`: Pass 2 on all intermediate drawings, deduplicated (`final.jsonl`,
  `final_<k>.json/.graphml` for the first 1000),
- and prints how many intermediate drawings can be finished.

## Monitoring

    squeue -u $USER -n nested_split
    cat runs/c12_run1/t*/done_units.txt | wc -l        # finished prefixes (compare with EXPECTED)
    grep PROGRESS logs/nested_split_<job>_<task>.log | tail -1

## Records of the 2026-10 runs
`JOBS.txt` (submitted jobs), `logs/` (Slurm logs) and `cluster_split_results_2026-10-05.tar.gz` are from the runs in
`/cluster/scratch/drivera/claude_split` (paths and script names there differ from the ones above; the final results
are the `c12_gk_l0` / `c14_gk_l0` runs of 2026-10-05).
