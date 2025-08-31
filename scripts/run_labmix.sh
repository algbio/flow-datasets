#!/usr/bin/env bash
set -euo pipefail

# Parameters
GENOMES=(1 2 3 4 5)
WINDOWS=(7000 10000)
K=15
DATASET="labmix"
EXTRA_OPTS="--mincycles 1 --nreads 5 --readlength 1000"
# Error epsilon variants: "none" means omit -e flag (uses default 1.0); others are passed to -e
E_VARIANTS=(none 0.75)
MAX_JOBS="${MAX_JOBS:-4}"  # override concurrently via: MAX_JOBS=8 ./scripts/run_ecoli_k63.sh
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT_DIR="${SCRIPT_DIR}/.."
cd "$ROOT_DIR"

# Create a FIFO-based semaphore with MAX_JOBS tokens
tmpfifo="$(mktemp -u)"
mkfifo "$tmpfifo"
exec 9<>"$tmpfifo"
rm -f "$tmpfifo"
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
        # perfect weights (no explicit -e provided)
        outdir="cyclic-graphs/${DATASET}/perfect-weights/g${g}-w${w}-k${K}/"
        e_flag=""
      else
        # imperfect weights: include epsilon value and place under imperfect-weights
        # replicate naming convention seen for k27 runs: append -e<value>
        outdir="cyclic-graphs/${DATASET}/imperfect-weights/g${g}-w${w}-k${K}-e${e}/"
        e_flag="-e ${e}"
      fi
      if [ -d "$outdir" ]; then
        echo "[skip] $outdir exists"
        continue
      fi
      read -u9  # acquire semaphore slot
      {
        cmd=(python3 construct.py -o "$outdir" -g "$g" -w "$w" -k "$K" $EXTRA_OPTS -D "$DATASET")
        if [ -n "$e_flag" ]; then
          cmd+=( $e_flag )
        fi
        echo "[run ] ${cmd[*]}"
        if "${cmd[@]}"; then
          echo "[done] $outdir"
        else
          echo "[fail] $outdir" >&2
          FAIL=1
        fi
        echo >&9  # release slot
      } &
      PIDS+=("$!")
    done
  done
done

wait  # wait for all background jobs

if [ $FAIL -ne 0 ]; then
  echo "One or more jobs failed." >&2
  exit 1
fi

echo "All jobs done."


