#!/bin/bash
# PROJECT: Wavefront Alignments Algorithms (Benchmarks)
# LICENCE: MIT License
# DESCRIPTION: A/B benchmark of the wavefront slab clear. A few long alignments
#              grow the slab, then many short alignments reuse the same aligner.
#              Each slab clear then sees a large slab and a small busy set.
# USAGE: ./wfa.bench.slab.sh BASE_BIN NEW_BIN [REPS]
#        BASE_BIN and NEW_BIN are align_benchmark binaries of the two builds.

# Config
PREFIX=$(dirname $0)
BASE=$1
NEW=$2
REPS=${3:-5}
GENERATE="$PREFIX/../bin/generate_dataset"

if [[ ! -x "$BASE" || ! -x "$NEW" ]]
then
    echo "USAGE: $0 BASE_BIN NEW_BIN [REPS]"
    exit -1
fi
if [ ! -x "$GENERATE" ]
then
    echo "[Error] $GENERATE not built. Please run make"
    exit -1
fi

# Generate the mixed dataset (10 x 10 kbp, then 200k x 100 bp)
DATA=$(mktemp -d)
trap 'rm -rf "$DATA"' EXIT
$GENERATE -o $DATA/long.seq -n 10 -l 10000 --error 0.05 > /dev/null
$GENERATE -o $DATA/short.seq -n 200000 -l 100 --error 0.05 > /dev/null
cat $DATA/long.seq $DATA/short.seq > $DATA/mixed.seq

# Print the alignment time of one run in ms
run() {
    "$1" -a gap-affine-wfa -i $DATA/mixed.seq --wfa-memory $2 2>&1 | \
        awk '/Time.Alignment/ {
            f = 1
            if ($4 == "m") f = 60000
            if ($4 == "s") f = 1000
            if ($4 == "us") f = 0.001
            if ($4 == "ns") f = 0.000001
            printf "%.2f\n", $3 * f
        }'
}

# Print the median of the arguments
median() {
    printf "%s\n" "$@" | sort -n | awk '{a[NR]=$1} END {print a[int((NR+1)/2)]}'
}

# Run both builds in turns, so that drift affects both the same way
printf "%-5s %12s %12s %8s\n" mem base_ms new_ms change
for mem in high med
do
    b=(); n=()
    for ((r=0;r<REPS;++r))
    do
        b+=($(run "$BASE" $mem))
        n+=($(run "$NEW" $mem))
    done
    mb=$(median "${b[@]}")
    mn=$(median "${n[@]}")
    change=$(awk -v b=$mb -v n=$mn 'BEGIN {printf "%+.1f%%", (n-b)/b*100}')
    printf "%-5s %12s %12s %8s\n" $mem $mb $mn $change
done
