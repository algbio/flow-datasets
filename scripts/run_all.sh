#!/usr/bin/env bash
# Run all dataset generation scripts sequentially with a unified MAX_JOBS setting.
# Usage:
#   ./scripts/run_all.sh            # uses default 32
#   MAX_JOBS=64 ./scripts/run_all.sh  # override
#
# Each underlying script already implements its own internal parallelism, so we
# keep them sequential here to avoid exploding total concurrent jobs.

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT_DIR="${SCRIPT_DIR}/.."
cd "$ROOT_DIR"

# Default if caller doesn't provide one
MAX_JOBS_VALUE="${MAX_JOBS:-32}"

echo "Using MAX_JOBS=${MAX_JOBS_VALUE} for all scripts"

# Collect run_*.sh scripts (excluding this file itself)
SCRIPTS=()
for f in "${SCRIPT_DIR}"/run_*.sh; do
  [ -e "$f" ] || continue
  base="$(basename "$f")"
  if [ "$base" = "run_all.sh" ]; then
    continue
  fi
  # Ensure executable bit
  [ -x "$f" ] || chmod +x "$f"
  SCRIPTS+=("$f")
done

if [ ${#SCRIPTS[@]} -eq 0 ]; then
  echo "No run_*.sh scripts found to execute." >&2
  exit 1
fi

echo "Found ${#SCRIPTS[@]} scripts:" >&2
for s in "${SCRIPTS[@]}"; do
  echo "  - $(basename "$s")" >&2
done

echo
FAIL=0
for s in "${SCRIPTS[@]}"; do
  echo "== Running $(basename "$s") =="
  # Export MAX_JOBS only for this invocation to avoid polluting environment after the loop.
  if MAX_JOBS="${MAX_JOBS_VALUE}" "$s"; then
    echo "[ok] $(basename "$s")"
  else
    echo "[fail] $(basename "$s")" >&2
    FAIL=1
  fi
  echo
done

if [ $FAIL -ne 0 ]; then
  echo "One or more scripts failed." >&2
  exit 1
fi

echo "All scripts completed successfully."
