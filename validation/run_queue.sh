#!/bin/bash
# usage: run_queue.sh <queue file> <log file> [timeout seconds per run, default 7000]
# each line of the queue: <binary> <args...>; runs sequentially, logs output and wall-clock time of every run
cd "$(dirname "$0")"
Q=$1; L=$2; TO=${3:-7000}
while read -r line; do
  [ -z "$line" ] && continue; [[ "$line" == \#* ]] && continue
  echo "=== $line   (start $(date '+%F %T'), timeout ${TO}s)" >> "$L"
  s=$(date +%s.%N)
  timeout "$TO" $line >> "$L" 2>&1; rc=$?
  e=$(date +%s.%N)
  echo "=== exit $rc wall $(echo "$e - $s" | bc) s" >> "$L"
done < "$Q"
echo "=== QUEUE DONE $(date '+%F %T')" >> "$L"
