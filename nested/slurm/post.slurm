#!/bin/bash -l
# After a Pass 1 array: merge its children, run Pass 2 on the base drawing and on all intermediate drawings,
# dedup the final drawings and count the intermediate drawings that can be finished.
# Submit from the repo root, after the array (job id J):
#   sbatch --dependency=afterany:J --export=ALL,BIN=.bin/split_c12_run1,RUN=c12_run1,EXPECTED=876629 nested/slurm/post.slurm
# Variables: BIN (same binary as the array), RUN, OUT (default runs), EXPECTED (number of prefixes, from
#   "BIN count DEPTH 1" with the same options; the script stops if Pass 1 is not complete)
#SBATCH --job-name=nested_post
#SBATCH --output=logs/%x_%j.log
#SBATCH --error=logs/%x_%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=8G
#SBATCH --time=24:00:00
#SBATCH --partition=standard
cd "$SLURM_SUBMIT_DIR"
module load stack/2024-06 boost >/dev/null 2>&1
set -e
OUT=${OUT:-runs}
R=$OUT/$RUN
DONE=$(cat $R/t*/done_units.txt | sort -u | wc -l)
echo "start $(date -Is), $RUN prefixes done: $DONE / $EXPECTED"
if [ "$DONE" -ne "$EXPECTED" ]; then echo "Pass 1 not complete: resubmit the $RUN array with the same settings, then this script."; exit 1; fi
"$BIN" merge ${R}_merged $R/t*                                                     # distinct intermediate drawings
"$BIN" search 0 1 0 2 ${R}_base_p2 > ${R}_base_p2.log                              # Pass 2 on the base drawing
"$BIN" merge_solutions ${R}_final_2cycles ${R}_base_p2
PROGRESS_SECONDS=600 "$BIN" search 0 1 0 2 ${R}_l1_p2 ${R}_merged/parents.txt > ${R}_l1_p2.log   # Pass 2 on all intermediate drawings
"$BIN" merge_solutions ${R}_final_3cycles ${R}_l1_p2
echo "finishable intermediate drawings: $(python3 -c "import json; print(len(set(json.loads(l)['meta']['parent'] for l in open('${R}_l1_p2/solutions.jsonl'))))" 2>/dev/null || echo 0)"
echo "end $(date -Is)"
