#!/bin/bash
# PROJECT: Wavefront Alignments Algorithms (Singletrack tests)
# LICENCE: MIT License
# AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
# DESCRIPTION: Focused Singletrack correctness tests
# USAGE: ./tests/wfa.singletrack.utest.sh [bin-dir]

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT_DIR="$(cd "$SCRIPT_DIR/.." && pwd)"
INPUT="$SCRIPT_DIR/wfa.utest.seq"
WORKDIR="$(mktemp -d "${TMPDIR:-/tmp}/wfa.singletrack.XXXXXX")"
ENDSFREE_INPUT="$WORKDIR/endsfree.seq"
EXTENSION_INPUT="$WORKDIR/extension.seq"

BIN="${1:-}/align_benchmark"
if [[ ! -f "$BIN" ]]; then
  BIN="$ROOT_DIR/bin/align_benchmark"
fi
if [[ ! -f "$BIN" ]]; then
  BIN="$ROOT_DIR/build/align_benchmark"
fi
if [[ ! -f "$BIN" ]]; then
  echo "[Error] Binaries not built. Please run cmake or make"
  exit 1
fi

cat > "$ENDSFREE_INPUT" <<'EOF'
>AAAACCCCGGGG
<TTTTAAAACCCCGGGGAAAA
>AAAACCCCGGGGTTTT
<GGGGAAAACCCCGGGG
>CCCCAAAAGGGG
<AAAAGGGGTTTT
EOF

cat > "$EXTENSION_INPUT" <<'EOF'
>AAAACCCCGGGGTTTT
<AAAACCCCGGGGAAAA
>ACGTACGTGGGG
<ACGTACGTTTTT
>TTTTAAAACCCC
<TTTTAAAAGGGG
EOF

compare_scores() {
  local base_out="$1"
  local singletrack_out="$2"
  awk '{print $1}' "$base_out" > "$base_out.score"
  awk '{print $1}' "$singletrack_out" > "$singletrack_out.score"
  diff -u "$base_out.score" "$singletrack_out.score"
}

run_case_input() {
  local name="$1"
  local input="$2"
  local algorithm="$3"
  shift 3

  local base_out="$WORKDIR/$name.high.alg"
  local singletrack_out="$WORKDIR/$name.singletrack.alg"

  echo ">>> Singletrack '$name'"
  "$BIN" -i "$input" -o "$base_out" -a "$algorithm" \
      --wfa-memory=high --check=correct "$@" > "$WORKDIR/$name.high.log" 2>&1
  "$BIN" -i "$input" -o "$singletrack_out" -a "$algorithm" \
      --wfa-memory=singletrack --check=correct "$@" > "$WORKDIR/$name.singletrack.log" 2>&1
  compare_scores "$base_out" "$singletrack_out"
}

run_case() {
  local name="$1"
  local algorithm="$2"
  shift 2

  run_case_input "$name" "$INPUT" "$algorithm" "$@"
}

run_span_heuristic_cases() {
  local span_name="$1"
  local input="$2"
  local span_arg="$3"

  run_case_input "affine.$span_name.wfadaptive" "$input" "gap-affine-wfa" \
      "$span_arg" "--wfa-heuristic=wfa-adaptive" "--wfa-heuristic-parameters=10,50,1"
  run_case_input "affine2p.$span_name.wfadaptive" "$input" "gap-affine2p-wfa" \
      "$span_arg" "--wfa-heuristic=wfa-adaptive" "--wfa-heuristic-parameters=10,50,1"
  run_case_input "affine.$span_name.xdrop" "$input" "gap-affine-wfa" \
      "$span_arg" "--wfa-heuristic=xdrop" "--wfa-heuristic-parameters=10,1"
  run_case_input "affine2p.$span_name.xdrop" "$input" "gap-affine2p-wfa" \
      "$span_arg" "--wfa-heuristic=xdrop" "--wfa-heuristic-parameters=10,1"
  run_case_input "affine.$span_name.zdrop" "$input" "gap-affine-wfa" \
      "$span_arg" "--wfa-heuristic=zdrop" "--wfa-heuristic-parameters=10,1"
  run_case_input "affine2p.$span_name.zdrop" "$input" "gap-affine2p-wfa" \
      "$span_arg" "--wfa-heuristic=zdrop" "--wfa-heuristic-parameters=10,1"
}

run_case "affine" "gap-affine-wfa"
run_case "affine2p" "gap-affine2p-wfa"
run_case "affine.p0" "gap-affine-wfa" "--affine-penalties=0,1,2,1"
run_case "affine.p1" "gap-affine-wfa" "--affine-penalties=0,3,1,4"
run_case "affine.p2" "gap-affine-wfa" "--affine-penalties=0,5,3,2"
run_case "affine.p3" "gap-affine-wfa" "--affine-penalties=-5,1,2,1"
run_case "affine.p4" "gap-affine-wfa" "--affine-penalties=-2,3,1,4"
run_case "affine.p5" "gap-affine-wfa" "--affine-penalties=-3,5,3,2"
run_case "affine.wfadaptive0" "gap-affine-wfa" \
    "--wfa-heuristic=wfa-adaptive" "--wfa-heuristic-parameters=10,50,1"
run_case "affine.wfadaptive1" "gap-affine-wfa" \
    "--wfa-heuristic=wfa-adaptive" "--wfa-heuristic-parameters=10,50,10"
run_case "affine.xdrop" "gap-affine-wfa" \
    "--wfa-heuristic=xdrop" "--wfa-heuristic-parameters=10,1"
run_case "affine.zdrop" "gap-affine-wfa" \
    "--wfa-heuristic=zdrop" "--wfa-heuristic-parameters=10,1"
run_case_input "affine.endsfree" "$ENDSFREE_INPUT" "gap-affine-wfa" "--wfa-span=ends-free,4,4,4,4"
run_case_input "affine2p.endsfree" "$ENDSFREE_INPUT" "gap-affine2p-wfa" "--wfa-span=ends-free,4,4,4,4"
run_case_input "affine.endsfree.negmatch" "$ENDSFREE_INPUT" "gap-affine-wfa" \
    "--wfa-span=ends-free,4,4,4,4" "--affine-penalties=-5,1,2,1"
run_case_input "affine2p.endsfree.negmatch" "$ENDSFREE_INPUT" "gap-affine2p-wfa" \
    "--wfa-span=ends-free,4,4,4,4" "--affine2p-penalties=-5,1,2,1,5,1"
run_span_heuristic_cases "endsfree" "$ENDSFREE_INPUT" "--wfa-span=ends-free,4,4,4,4"
run_case_input "affine.extension" "$EXTENSION_INPUT" "gap-affine-wfa" "--wfa-span=extension"
run_case_input "affine2p.extension" "$EXTENSION_INPUT" "gap-affine2p-wfa" "--wfa-span=extension"
run_case_input "affine.extension.negmatch" "$EXTENSION_INPUT" "gap-affine-wfa" \
    "--wfa-span=extension" "--affine-penalties=-5,1,2,1"
run_case_input "affine2p.extension.negmatch" "$EXTENSION_INPUT" "gap-affine2p-wfa" \
    "--wfa-span=extension" "--affine2p-penalties=-5,1,2,1,5,1"
run_span_heuristic_cases "extension" "$EXTENSION_INPUT" "--wfa-span=extension"

"$SCRIPT_DIR/wfa.singletrack.proof.py" "$BIN"

expect_reject() {
  local name="$1"
  shift
  if "$BIN" -i "$INPUT" -o "$WORKDIR/$name.alg" "$@" > "$WORKDIR/$name.log" 2>&1; then
    echo "[Error] Singletrack accepted unsupported configuration '$name'"
    exit 1
  fi
}

expect_reject "score-only" -a gap-affine-wfa --wfa-memory=singletrack --wfa-score-only
expect_reject "edit" -a edit-wfa --wfa-memory=singletrack
expect_reject "lambda" -a gap-affine-wfa --wfa-memory=singletrack --wfa-lambda
expect_reject "banded-static" -a gap-affine-wfa --wfa-memory=singletrack \
    --wfa-heuristic=banded-static --wfa-heuristic-parameters=-10,10
expect_reject "banded-adaptive" -a gap-affine-wfa --wfa-memory=singletrack \
    --wfa-heuristic=banded-adaptive --wfa-heuristic-parameters=-10,10,1

echo ">>> Singletrack tests passed ($WORKDIR)"
