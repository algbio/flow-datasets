#!/usr/bin/env bash
set -euo pipefail

# Parameters (adjust as desired)
GENOMES=(5 10 15 20 25 30 40 50)  # cap at downloaded genomes (max 100)
WINDOWS=(10000)  # SARS-CoV-2 ~29.9kb
K=27
DATASET="sars_cov_2"
EXTRA_OPTS="--mincycles 1 --nreads 5 --readlength 1500"
E_VARIANTS=(none imp)
MAX_JOBS="${MAX_JOBS:-6}"
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT_DIR="${SCRIPT_DIR}/.."
cd "$ROOT_DIR"

# simple job slot semaphore
fifo="$(mktemp -u)"
mkfifo "$fifo"
exec 9<>"$fifo"
rm -f "$fifo"
for ((i=0;i<MAX_JOBS;i++)); do echo >&9; done

FAIL=0
PIDS=()

cleanup() {
  trap - INT TERM EXIT
  if [ ${#PIDS[@]} -gt 0 ]; then
    kill ${PIDS[@]} 2>/dev/null || true
  fi
  exec 9>&- 9<&- || true
}
trap cleanup INT TERM EXIT

for g in "${GENOMES[@]}"; do
  for w in "${WINDOWS[@]}"; do
    for e in "${E_VARIANTS[@]}"; do
      if [ "$e" = none ]; then
        outdir="cyclic-graphs/${DATASET}/perfect-weights/g${g}-w${w}-k${K}/"
        e_flag=""
      else
        outdir="cyclic-graphs/${DATASET}/imperfect-weights/g${g}-w${w}-k${K}-imp/"
        e_flag="--poisson-edge-errors"
      fi
      if [ -d "$outdir" ]; then
        echo "[skip] $outdir exists"
        continue
      fi
      read -u9
      {
        cmd=(python3 construct.py -o "$outdir" -g "$g" -w "$w" -k "$K" $EXTRA_OPTS -D "$DATASET")
        if [ -n "$e_flag" ]; then
          cmd+=( $e_flag )
        fi
        echo "[run ] ${cmd[*]}"
        if ${cmd[@]}; then
          echo "[done] $outdir"
        else
          echo "[fail] $outdir" >&2
          FAIL=1
        fi
        echo >&9
      } &
      PIDS+=("$!")
    done
  done
  # Don't request more genomes than available FASTAs
  # (User should manually trim GENOMES list if fewer downloaded.)
  # We only warn after first iteration if g exceeds available.
  avail=$(ls datasets/${DATASET}/*.fasta.gz 2>/dev/null | wc -l | tr -d ' ' || echo 0)
  if [ "$avail" -gt 0 ] && [ "$g" -gt "$avail" ]; then
    echo "[warn] Requested g=$g exceeds available $avail FASTA files; further jobs may fail." >&2
  fi

done

wait

if [ $FAIL -ne 0 ]; then
  echo "One or more jobs failed." >&2
  exit 1
fi

echo "All jobs done."
