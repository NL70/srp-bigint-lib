#!/bin/bash
set -e

cd /app
mkdir -p results

CPU_CORE=0
TIMESTAMP=$(date +%Y%m%d_%H%M%S)
LOGFILE="results/run_${TIMESTAMP}.log"
CSVFILE="results/run_${TIMESTAMP}.csv"
REPS=50

echo "Running testBigInt pinned to CPU core $CPU_CORE with $REPS repetitions..."

# Run the benchmark and save full output to log
taskset -c "$CPU_CORE" ./bin/testBigInt \
  --benchmark_repetitions="$REPS" \
  --benchmark_enable_random_interleaving=true \
  --benchmark_format=csv \
  2>&1 | tee "$LOGFILE"

# Extract only the CSV lines (keep header info in .log)
tail -n +5 "$LOGFILE" > "$CSVFILE"
echo "Benchmark completed."
echo "Full log saved to: $LOGFILE"
echo "CSV data saved to: $CSVFILE"
