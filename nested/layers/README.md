# All drawings of the nested cycle graphs G_{12,k} and G_{14,k}

G_{c,k}: cycles D1, ..., Dk of length c (12 or 14), consecutive cycles connected by the braided matching
(even i -> (i-3)', odd i -> (i+3)'). D1 and Dk are uncrossed (in the full graph the gadgets force this),
D2, ..., D(k-1) may be crossed. The code finds all 3-plane drawings of G_{c,k} and tests them for saturation.

## Files

- `nested_c12_layers.cpp`, `nested_c14_layers.cpp`: the setup, everything hard-coded at the top of `main`: the cycle
  length, the number of cycles k, the output directory `run_dir`, the number of tasks and the split depth, the layer
  edges (new cycle + matching, written out), and the checkpoints of the flow check. Change `num_cycles` for
  more cycles.
- `nested_layers.h`: the driver. Builds the drawings cycle by cycle with the searcher, splits the work into tasks,
  tests the complete drawings and combines the tasks.
- `nested_search.h`: the searcher. Adds **one** layer (the next cycle and its matching) to a given drawing, in all
  possible ways. Includes all the heuristic code.
- `layers.slurm`: one task as a Slurm array task.

## The searcher (nested_search.h)

`LayerSearcher` gets a drawing whose last ("active") cycle is uncrossed or crossable, and adds the next layer by a
depth-first search over the possible paths (at most 3 crossings) of the layer edges, one edge after the other:

- **Pass 1**: the new cycle is crossable. Every complete drawing that still has a *potential final face* (a face in
  which a further uncrossed cycle can be drawn, i.e. flow >= c from the new cycle vertices) is an *intermediate drawing* (a leaf of the
  search tree). It is reduced to its relevant part (as in Hoffmann et al.: the faces the next layer can use, the rest replaced by
  uncrossable stars) and reported together with the canonical code of the reduced drawing.
- **Pass 2**: the new cycle is uncrossed. Every complete drawing is reported.

Pruning:

- reachability: every remaining edge must still have a path on its own; the remaining edges are reordered, fewest
  options first;
- multicommodity flow (at the checkpoints): the remaining edges must be routable together;
- Pass 1 only: a potential final face must still be possible, every unplaced new vertex must still have a face it
  can go in, and paths that put a new vertex into an impossible face are skipped.

Canonical codes: two drawings get the same code iff they are isomorphic (planarizations, possibly mirrored, with
crossing counts). `canonical_code` fixes the active cycle up to the symmetries of the layer edges (used for the
intermediate drawings, which are extended further by label), `canonical_code_free` has no fixed labels (used for
the complete drawings).

Splitting: the paths of the first `split_depth` edges ("units") are numbered in DFS order, identically in every
task; a task only searches below the units with `index % num_tasks == task_index`. Finished units are written to
`done_units.txt` and skipped when a killed task is restarted.

## The driver (nested_layers.h)

`run <task_index>` (one task):

1. Pass 1 on the base drawing (D1 with an uncrossable star inside) adds the crossable D2, only below this task's
   units.
2. At every intermediate drawing with j cycles (distinct within the task): Pass 1 again for D(j+1) as long as j+1 < k, Pass 2 for
   the uncrossed Dk otherwise. With k = 2 there is no Pass 1: Pass 2 adds D2 to D1 directly (task 0 only).
3. Every complete drawing (distinct within the task) is tested for saturation: an uncrossable star is added in the
   outer face of Dk, and every missing edge between two cycles, or within a middle cycle, is tried with all routes
   of at most 3 crossings (edges within D1 and within Dk count as present). The drawing is written to
   `<run_dir>/t<task_index>/c12_3cycles_complete.jsonl` (etc.) with its insertable edges.
4. The hashes of the intermediate drawings (and of those with k-1 cycles that could be finished) go to
   `intermediate.txt`, needed to count distinct drawings over all tasks.

`report` (after all tasks): combines the tasks and prints

- the distinct intermediate drawings with j cycles (the last one crossable), j = 2, ..., k-1,
- how many intermediate drawings with k-1 cycles can be finished with an uncrossed Dk,
- the distinct complete drawings, how many are saturated, and which kinds of edges are insertable,

and writes all distinct complete drawings to `<run_dir>/c12_3cycles_complete.jsonl`.

Results with k = 3 (thesis): C12: 16,611 intermediate drawings with two cycles, 1,517 of them finishable, 3,231
complete drawings; C14: 44,682, 6,311 and 15,445.

## Running

From the repo root (`run_dir` is relative to the directory the binary is run from):

    mkdir -p .bin logs
    g++ -O3 -std=c++17 -Iinclude nested/layers/nested_c12_layers.cpp -o .bin/nested_c12_layers
    sbatch --array=0-99 --time=04:00:00 --export=ALL,BIN=.bin/nested_c12_layers nested/layers/layers.slurm
    .bin/nested_c12_layers report     # after the array

or locally, e.g. `seq 0 99 | xargs -P 8 -n 1 .bin/nested_c12_layers run`. The array size must equal `num_tasks`.
Killed tasks: rerun the same command (finished units are skipped).
With k = 3 and 100 tasks: C12 takes about 5-10 min per task, C14 about four times as long.
