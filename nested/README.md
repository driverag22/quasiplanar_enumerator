# Nested cycles C12 / C14

- `layers/`: the thesis code. All drawings of G_{c,k} for c = 12, 14 (k cycles, the first and last uncrossed) and
  the counts of the thesis (16,611 / 44,682 drawings with two cycles, the second crossable, etc.).
  `nested_c12_layers.cpp` / `nested_c14_layers.cpp` (setup, everything hard-coded), `nested_layers.h` (driver),
  `nested_search.h` (search for one layer).
- `exploration/`: the code of the explorations (extra edge classes, open runs, sweeps) and the original split
  pipeline (merge, Pass 2 on given parent drawings, e.g. the 17-parent sample with four cycles).
  `iterative_nested_c12_split.cpp` / `_c14_split.cpp` (gadget and extra classes), `iterative_nested_split.h`
  (search, with options for extra edges), `split_driver.h` (command line: `search`, `count`, `merge`,
  `merge_solutions`, `saturate`), Slurm scripts in `exploration/slurm/`.
- `oldCode/`: the original (serial) code.
- `smallTests/`: small single-drawing tests.

All commands below are run from the **repo root**; the Slurm scripts `cd` to the directory `sbatch` was called
from. Slurm writes its logs to `logs/` there, which has to exist: `mkdir -p logs`. On Euler, build on the login node
after `bash -l; module load stack/2024-06 boost`.

## Thesis code (layers/)
The number of cycles, the output directory `run_dir`, the number of tasks and the split depth are fixed at the top
of `main` in `nested_c12_layers.cpp` / `nested_c14_layers.cpp`.

    mkdir -p .bin
    g++ -O3 -std=c++17 -Iinclude nested/layers/nested_c12_layers.cpp -o .bin/nested_c12_layers
    sbatch --array=0-99 --time=04:00:00 --export=ALL,BIN=.bin/nested_c12_layers nested/layers/layers.slurm
    .bin/nested_c12_layers report     # after the array

Without a cluster the tasks can run in parallel locally, e.g. `seq 0 99 | xargs -P 8 -n 1 .bin/nested_c12_layers run`.
The report prints the distinct intermediate drawings per number of cycles, how many intermediate drawings with
k-1 cycles can be finished and the distinct complete drawings, and writes all distinct
complete drawings to `<run_dir>/c12_3cycles_complete.jsonl` (etc.). Killed tasks: rerun the same command
(finished prefixes are skipped).

## Exploration code: split pipeline (exploration/)

### 1. Build
On the login node (or `sbatch nested/exploration/slurm/build.slurm`):

    mkdir -p .bin
    g++ -O3 -Wall -std=c++17 -Iinclude nested/exploration/iterative_nested_c12_split.cpp -o .bin/split_c12

Freeze a copy per run (`cp .bin/split_c12 .bin/split_c12_run1`) and use it for all tasks, resubmissions and the
post job of that run: the prefix numbering must be identical in every task.

### 2. Choose the split depth
The options are read from the environment; set them exactly as for the run:

    export REACH_EVERY_EDGE_FROM=1 FINAL_FACE_EVERY_EDGE_FROM=1 PLACEMENT_EVERY_EDGE_FROM=0 PLACEMENT_FILTER=1
    .bin/split_c12_run1 count 12 1        # prints the number of prefixes at depth 12 (= EXPECTED below)

With these options: C12 depth 12 -> 876,629 prefixes, C14 depth 12 -> 877,520.

### 3. Pass 1 on the base drawing (job array)

    sbatch --array=0-99 --time=04:00:00 \
      --export=ALL,BIN=.bin/split_c12_run1,RUN=c12_run1,TASKS=100,DEPTH=12,PASSES=1,REACH_EVERY_EDGE_FROM=1,FINAL_FACE_EVERY_EDGE_FROM=1,PLACEMENT_EVERY_EDGE_FROM=0,PLACEMENT_FILTER=1,PROGRESS_SECONDS=1800 \
      nested/exploration/slurm/search.slurm

Keep the requests close to the real usage (level 0: at most ~35 min and ~155 MB per task): small jobs start
sooner, and the per-user memory limit (128 GB) otherwise caps the number of tasks running at once.

Task i searches below the prefixes j with j mod TASKS == i and writes `runs/c12_run1/t<i>/`:
`children.jsonl` (intermediate drawings, deduplicated per task only), `done_units.txt` (finished prefixes),
`extensible.txt`. If tasks are killed (time limit), resubmit the same command: finished prefixes are skipped.

### 4. Merge + Pass 2 (one job, after the array)

    sbatch --dependency=afterany:<array job id> \
      --export=ALL,BIN=.bin/split_c12_run1,RUN=c12_run1,EXPECTED=876629 nested/exploration/slurm/post.slurm

Checks that all prefixes are done, then writes
- `runs/c12_run1_merged/parents.jsonl`: distinct intermediate drawings (gadget-aware isomorphism),
- `runs/c12_run1_final_2cycles/`: Pass 2 on the base drawing, deduplicated (strong isomorphism),
- `runs/c12_run1_final_3cycles/`: Pass 2 on all intermediate drawings, deduplicated (`final.jsonl`,
  `final_<k>.json/.graphml` for the first 1000),
- and prints how many intermediate drawings can be finished.

### Sweeps of extra edge classes (SAT_LAYERS / SAT_CHORDS, split_driver.h)
`exploration/slurm/layers_sweep.slurm`: one array over several variants, index = variant * TASKS + task. The
variants file `<scratch>/runs/<SWEEP>/variants.txt` has one line `<name> <SAT_LAYERS> <SAT_CHORDS>` per variant;
`FAILFAST=1` stops a variant at its first drawing with an insertable edge (marker `<name>/not_maximal`).

    sbatch --array=0-863%10 --export=ALL,BIN=split_c12_v17,SWEEP=screw4ff,TASKS=16,FAILFAST=1 layers_sweep.slurm

### Monitoring

    squeue -u $USER -n nested_split
    cat runs/c12_run1/t*/done_units.txt | wc -l        # finished prefixes (compare with EXPECTED)
    grep PROGRESS logs/nested_split_<job>_<task>.log | tail -1

### Records of the 2026-10 runs
`exploration/slurm/JOBS.txt` (submitted jobs), `exploration/slurm/logs/` (Slurm logs) and
`cluster_split_results_2026-10-05.tar.gz` are from the runs in `/cluster/scratch/drivera/claude_split` (paths and
script names there differ from the ones above; the final results are the `c12_gk_l0` / `c14_gk_l0` runs of
2026-10-05).
