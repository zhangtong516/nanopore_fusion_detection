#!/usr/bin/env bash
set -euo pipefail
sample=${1:?}
shift
if [ $# -eq 0 ]; then echo "No FASTQ inputs" >&2; exit 1; fi
tmp_list=()
for f in "$@"; do
  tmp_list+=("$f")
done
for f in "${tmp_list[@]}"; do
  if [[ "$f" =~ \.gz$ ]]; then zcat "$f"; else cat "$f"; fi
done | gzip -c > ${sample}_merged.fastq.gz

