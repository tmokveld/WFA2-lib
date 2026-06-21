#!/usr/bin/env python3
"""BiWFA ends-free zero-match proof checks."""

from __future__ import annotations

from dataclasses import dataclass
from itertools import product
import random
import re
import subprocess
import sys
import tempfile
from pathlib import Path


DNA = "ACGT"
TINY_DNA = "AC"
TINY_MAX_LENGTH = 3
CIGAR_RE = re.compile(r"(\d*)([=MXID])")
NEGATIVE_MATCH_MESSAGE = (
    "[WFA] BiWFA ends-free with negative match rewards (match < 0) "
    "requires aligned-length-aware breakpoint scoring and is not implemented yet\n"
)


@dataclass(frozen=True)
class Case:
    pattern: str
    text: str
    terminal: tuple[str, str, int] | None = None


@dataclass(frozen=True)
class AlignmentRow:
    score: int
    pattern: str
    text: str
    cigar: str


@dataclass(frozen=True)
class RunResult:
    rows: list[AlignmentRow]
    command: list[str]
    input_path: Path
    output_path: Path


METRICS = [
    ("indel", "indel-wfa", {"kind": "indel"}),
    ("edit", "edit-wfa", {"kind": "edit"}),
    (
        "linear",
        "gap-linear-wfa",
        {"kind": "linear", "mismatch": 4, "indel": 2},
    ),
    (
        "affine",
        "gap-affine-wfa",
        {"kind": "affine", "mismatch": 4, "open": 6, "extend": 2},
    ),
    (
        "affine2p",
        "gap-affine2p-wfa",
        {
            "kind": "affine2p",
            "mismatch": 4,
            "open1": 6,
            "extend1": 2,
            "open2": 24,
            "extend2": 1,
        },
    ),
]


FREE_END_CONFIGS = [
    ("pattern-begin", (4, 0, 0, 0)),
    ("pattern-end", (0, 4, 0, 0)),
    ("text-begin", (0, 0, 4, 0)),
    ("text-end", (0, 0, 0, 4)),
    ("pattern-only", (4, 4, 0, 0)),
    ("text-only", (0, 0, 4, 4)),
    ("overlap", (4, 0, 0, 4)),
    ("all", (4, 4, 4, 4)),
]


TINY_FREE_END_CONFIGS = [
    ("tiny-global", (0, 0, 0, 0)),
    # align_benchmark treats values in [0.0,1.0] as proportions, not counts.
    ("tiny-pattern-begin", (2, 0, 0, 0)),
    ("tiny-pattern-end", (0, 2, 0, 0)),
    ("tiny-text-begin", (0, 0, 2, 0)),
    ("tiny-text-end", (0, 0, 0, 2)),
    ("tiny-all", (2, 2, 2, 2)),
]


HEURISTIC_SMOKE_CONFIGS = [
    (
        "wfa-adaptive",
        ["--wfa-heuristic=wfa-adaptive", "--wfa-heuristic-parameters=10,50,1"],
    ),
    ("wfmash", ["--wfa-heuristic=wfmash", "--wfa-heuristic-parameters=10,50,1"]),
]


GLOBAL_CONFIG = ("global", (0, 0, 0, 0))


def find_align_benchmark(argv: list[str]) -> Path:
    if len(argv) > 1:
        candidate = Path(argv[1])
        if candidate.is_dir():
            candidate = candidate / "align_benchmark"
        if candidate.is_file():
            return candidate.resolve()
    root = Path(__file__).resolve().parent.parent
    for candidate in (root / "bin" / "align_benchmark", root / "build" / "align_benchmark"):
        if candidate.is_file():
            return candidate.resolve()
    raise SystemExit("[Error] Binaries not built. Please run cmake or make")


def random_dna(rng: random.Random, length: int) -> str:
    return "".join(rng.choice(DNA) for _ in range(length))


def mutate(rng: random.Random, sequence: str, edits: int) -> str:
    bases = list(sequence)
    for _ in range(edits):
        op = rng.choice("XXIID") if bases else "I"
        if op == "X":
            idx = rng.randrange(len(bases))
            bases[idx] = rng.choice([base for base in DNA if base != bases[idx]])
        elif op == "I":
            bases.insert(rng.randrange(len(bases) + 1), rng.choice(DNA))
        elif op == "D":
            del bases[rng.randrange(len(bases))]
    return "".join(bases) or rng.choice(DNA)


def divergent_pair(length: int) -> tuple[str, str]:
    rng = random.Random(0xB1FAD)
    pattern = random_dna(rng, length)
    text = list(pattern)
    for idx in range(0, len(text), 3):
        text[idx] = rng.choice([base for base in DNA if base != text[idx]])
    return pattern, "".join(text)


def enumerate_sequences(alphabet: str, max_length: int) -> list[str]:
    sequences: list[str] = []
    for length in range(1, max_length + 1):
        sequences.extend("".join(chars) for chars in product(alphabet, repeat=length))
    return sequences


def exhaustive_tiny_cases(free_ends: tuple[int, int, int, int]) -> list[Case]:
    sequences = enumerate_sequences(TINY_DNA, TINY_MAX_LENGTH)
    p0, pf, t0, tf = free_ends
    cases: list[Case] = []
    for pattern in sequences:
        if p0 > len(pattern) or pf > len(pattern):
            continue
        for text in sequences:
            if t0 > len(text) or tf > len(text):
                continue
            cases.append(Case(pattern, text))
    return cases


def generated_cases(free_ends: tuple[int, int, int, int]) -> list[Case]:
    p0, pf, t0, tf = free_ends
    rng = random.Random(0x5EED + p0 * 101 + pf * 17 + t0 * 7 + tf)
    core = "ACGTACGTGACCTGACTGA"
    cases = [
        Case("ACGTACGTACGT", "ACGTTCGTACGA"),
        Case("GGGGAAAACCCCTTTT", "GGGGAAATCCCCTTTT"),
        Case("TGCATGCATGCA", "TGCAGGCATGCA"),
    ]
    if p0:
        cases.append(Case("T" * p0 + core, core, ("prefix", "D", p0)))
        cases.append(Case("T" * (p0 + 2) + core, core))
    if pf:
        cases.append(Case(core + "G" * pf, core, ("suffix", "D", pf)))
        cases.append(Case(core + "G" * (pf + 2), core))
    if t0:
        cases.append(Case(core, "T" * t0 + core, ("prefix", "I", t0)))
        cases.append(Case(core, "T" * (t0 + 2) + core))
    if tf:
        cases.append(Case(core, core + "G" * tf, ("suffix", "I", tf)))
        cases.append(Case(core, core + "G" * (tf + 2)))
    if p0 and tf:
        cases.append(Case("C" * p0 + core, core + "A" * tf))
    if t0 and pf:
        cases.append(Case(core + "C" * pf, "G" * t0 + core))
    for _ in range(8):
        base = random_dna(rng, rng.randint(16, 38))
        mate = mutate(rng, base, rng.randint(1, 8))
        pattern = (
            random_dna(rng, rng.randint(0, p0))
            + base
            + random_dna(rng, rng.randint(0, pf))
        )
        text = (
            random_dna(rng, rng.randint(0, t0))
            + mate
            + random_dna(rng, rng.randint(0, tf))
        )
        cases.append(Case(pattern, text))
    long_pattern, long_text = divergent_pair(900)
    cases.append(
        Case(
            random_dna(rng, p0) + long_pattern + random_dna(rng, pf),
            random_dna(rng, t0) + long_text + random_dna(rng, tf),
        )
    )
    return cases


def write_cases(path: Path, cases: list[Case]) -> None:
    with path.open("w") as handle:
        for case in cases:
            print(f">{case.pattern}", file=handle)
            print(f"<{case.text}", file=handle)


def run_align_benchmark(
    align_benchmark: Path,
    workdir: Path,
    name: str,
    cases: list[Case],
    algorithm: str,
    memory_mode: str,
    span: str,
    score_only: bool,
    extra_args: list[str] | None = None,
) -> RunResult:
    input_path = workdir / f"{name}.seq"
    output_path = workdir / f"{name}.{memory_mode}.alg"
    write_cases(input_path, cases)
    cmd = [
        str(align_benchmark),
        "--input",
        str(input_path),
        "--output-full",
        str(output_path),
        "--algorithm",
        algorithm,
        f"--wfa-memory={memory_mode}",
        f"--wfa-span={span}",
    ]
    if score_only:
        cmd.append("--wfa-score-only")
    if extra_args:
        cmd.extend(extra_args)
    completed = subprocess.run(cmd, text=True, capture_output=True)
    if completed.returncode != 0:
        raise RuntimeError(
            "align_benchmark failed:\n"
            + " ".join(cmd)
            + "\n\n"
            + completed.stdout
            + completed.stderr
        )
    rows: list[AlignmentRow] = []
    with output_path.open() as handle:
        for line in handle:
            cols = line.rstrip("\n").split("\t")
            if len(cols) != 6:
                raise AssertionError(f"unexpected output row: {line!r}")
            _, _, score, pattern, text, cigar = cols
            rows.append(AlignmentRow(int(score), pattern, text, cigar))
    if len(rows) != len(cases):
        raise AssertionError(f"expected {len(cases)} outputs, got {len(rows)}")
    return RunResult(rows, cmd, input_path, output_path)


def parse_cigar(cigar: str) -> list[str]:
    if cigar == "-":
        return []
    ops: list[str] = []
    pos = 0
    for match in CIGAR_RE.finditer(cigar):
        if match.start() != pos:
            raise AssertionError(f"malformed CIGAR {cigar!r}")
        count = int(match.group(1) or "1")
        ops.extend(match.group(2) for _ in range(count))
        pos = match.end()
    if pos != len(cigar):
        raise AssertionError(f"malformed CIGAR {cigar!r}")
    return ops


def validate_cigar(pattern: str, text: str, cigar: str) -> list[tuple[int, int]]:
    v = 0
    h = 0
    coords = [(0, 0)]
    for op in parse_cigar(cigar):
        if op in ("M", "="):
            if v >= len(pattern) or h >= len(text) or pattern[v] != text[h]:
                raise AssertionError(f"invalid match at ({v},{h}) in {cigar}")
            v += 1
            h += 1
        elif op == "X":
            if v >= len(pattern) or h >= len(text) or pattern[v] == text[h]:
                raise AssertionError(f"invalid mismatch at ({v},{h}) in {cigar}")
            v += 1
            h += 1
        elif op == "I":
            if h >= len(text):
                raise AssertionError(f"insertion overruns text in {cigar}")
            h += 1
        elif op == "D":
            if v >= len(pattern):
                raise AssertionError(f"deletion overruns pattern in {cigar}")
            v += 1
        coords.append((v, h))
    if v != len(pattern) or h != len(text):
        raise AssertionError(
            f"CIGAR ends at ({v},{h}), expected ({len(pattern)},{len(text)}) in {cigar}"
        )
    return coords


def score_ops(ops: list[str], penalties: dict[str, int]) -> int:
    kind = penalties["kind"]
    if kind == "indel":
        return sum(1 for op in ops if op in ("I", "D")) + 2 * sum(1 for op in ops if op == "X")
    if kind == "edit":
        return sum(1 for op in ops if op in ("X", "I", "D"))
    if kind == "linear":
        score = 0
        for op in ops:
            if op == "X":
                score -= penalties["mismatch"]
            elif op in ("I", "D"):
                score -= penalties["indel"]
        return score
    score = 0
    idx = 0
    while idx < len(ops):
        op = ops[idx]
        end = idx + 1
        while end < len(ops) and ops[end] == op:
            end += 1
        length = end - idx
        if op == "X":
            score -= penalties["mismatch"] * length
        elif op in ("I", "D"):
            if kind == "affine":
                score -= penalties["open"] + penalties["extend"] * length
            elif kind == "affine2p":
                gap1 = penalties["open1"] + penalties["extend1"] * length
                gap2 = penalties["open2"] + penalties["extend2"] * length
                score -= min(gap1, gap2)
        idx = end
    return score


def cigar_endsfree_score(
    pattern: str,
    text: str,
    cigar: str,
    penalties: dict[str, int],
    free_ends: tuple[int, int, int, int],
) -> int:
    ops = parse_cigar(cigar)
    coords = validate_cigar(pattern, text, cigar)
    p0, pf, t0, tf = free_ends
    plen = len(pattern)
    tlen = len(text)
    starts = [
        idx
        for idx, (v, h) in enumerate(coords)
        if (h == 0 and v <= p0) or (v == 0 and h <= t0)
    ]
    ends = [
        idx
        for idx, (v, h) in enumerate(coords)
        if (v == plen and h >= tlen - tf) or (h == tlen and v >= plen - pf)
    ]
    candidates = [
        score_ops(ops[start:end], penalties)
        for start in starts
        for end in ends
        if end >= start
    ]
    if not candidates:
        raise AssertionError(f"CIGAR has no legal ends-free segment: {cigar}")
    if penalties["kind"] in ("indel", "edit"):
        return min(candidates)
    return max(candidates)


def fail_with_replay(
    message: str,
    high: RunResult,
    biwfa: RunResult,
    case_index: int,
    observed: int,
    expected: int,
) -> None:
    raise AssertionError(
        f"{message}\n"
        f"  case: {case_index}\n"
        f"  observed: {observed}\n"
        f"  expected: {expected}\n"
        f"  pattern: {high.rows[case_index].pattern}\n"
        f"  text: {high.rows[case_index].text}\n"
        f"  biwfa_cigar: {biwfa.rows[case_index].cigar}\n"
        f"  replay_high: {' '.join(high.command)}\n"
        f"  replay_biwfa: {' '.join(biwfa.command)}"
    )


def check_cases(
    align_benchmark: Path,
    workdir: Path,
    label: str,
    free_ends: tuple[int, int, int, int],
    span: str,
    cases: list[Case],
    metrics: list[tuple[str, str, dict[str, int]]] | None = None,
    extra_args: list[str] | None = None,
) -> None:
    metrics = METRICS if metrics is None else metrics
    for metric_label, algorithm, penalties in metrics:
        print(f">>> Proving BiWFA {label} {metric_label}")
        name = f"{label}.{metric_label}"
        high_score = run_align_benchmark(
            align_benchmark, workdir, f"{name}.high.score",
            cases, algorithm, "high", span, True, extra_args)
        biwfa_score = run_align_benchmark(
            align_benchmark, workdir, f"{name}.biwfa.score",
            cases, algorithm, "ultralow", span, True, extra_args)
        biwfa_full = run_align_benchmark(
            align_benchmark, workdir, f"{name}.biwfa.full",
            cases, algorithm, "ultralow", span, False, extra_args)
        for idx, (expected_row, score_row, full_row) in enumerate(
            zip(high_score.rows, biwfa_score.rows, biwfa_full.rows)
        ):
            if expected_row.pattern != score_row.pattern or expected_row.text != score_row.text:
                raise AssertionError(f"score-only output sequence changed in case {idx}")
            if expected_row.score != score_row.score:
                fail_with_replay(
                    "BiWFA score-only score differs from high-memory score-only",
                    high_score, biwfa_score, idx, score_row.score, expected_row.score)
            if expected_row.pattern != full_row.pattern or expected_row.text != full_row.text:
                raise AssertionError(f"full output sequence changed in case {idx}")
            observed = cigar_endsfree_score(
                full_row.pattern, full_row.text, full_row.cigar, penalties, free_ends)
            if observed != expected_row.score:
                fail_with_replay(
                    "BiWFA full CIGAR does not achieve high-memory ends-free score",
                    high_score, biwfa_full, idx, observed, expected_row.score)
            check_terminal_convention(cases[idx], full_row)


def check_terminal_convention(case: Case, row: AlignmentRow) -> None:
    if case.terminal is None:
        return
    side, op, count = case.terminal
    ops = parse_cigar(row.cigar)
    observed = ops[:count] if side == "prefix" else ops[-count:]
    expected = [op] * count
    if observed != expected:
        raise AssertionError(
            f"terminal free overhang not represented as {op}: {row.cigar}"
        )


def check_config(
    align_benchmark: Path,
    workdir: Path,
    label: str,
    free_ends: tuple[int, int, int, int],
    span: str,
) -> None:
    cases = generated_cases(free_ends)
    check_cases(align_benchmark, workdir, label, free_ends, span, cases)


def check_tiny_exhaustive(align_benchmark: Path, workdir: Path) -> None:
    for label, free_ends in TINY_FREE_END_CONFIGS:
        span = "global" if free_ends == (0, 0, 0, 0) else "ends-free,%d,%d,%d,%d" % free_ends
        check_cases(
            align_benchmark, workdir, label, free_ends, span,
            exhaustive_tiny_cases(free_ends))


def check_negative_match_rejection(align_benchmark: Path, workdir: Path) -> None:
    configs = [
        (
            "affine-score-only",
            "gap-affine-wfa",
            ["--affine-penalties=-5,4,6,2"],
            True,
        ),
        (
            "affine-full",
            "gap-affine-wfa",
            ["--affine-penalties=-5,4,6,2"],
            False,
        ),
        (
            "affine2p-score-only",
            "gap-affine2p-wfa",
            ["--affine2p-penalties=-5,4,6,2,24,1"],
            True,
        ),
        (
            "affine2p-full",
            "gap-affine2p-wfa",
            ["--affine2p-penalties=-5,4,6,2,24,1"],
            False,
        ),
    ]
    cases = [Case("ACGTACGT", "ACGTTACGT")]
    span = "ends-free,2,0,0,0"
    for label, algorithm, penalty_args, score_only in configs:
        print(f">>> Proving BiWFA rejects negative-match {label}")
        input_path = workdir / f"reject.{label}.seq"
        output_path = workdir / f"reject.{label}.alg"
        write_cases(input_path, cases)
        cmd = [
            str(align_benchmark),
            "--input",
            str(input_path),
            "--output-full",
            str(output_path),
            "--algorithm",
            algorithm,
            "--wfa-memory=ultralow",
            f"--wfa-span={span}",
            *penalty_args,
        ]
        if score_only:
            cmd.append("--wfa-score-only")
        completed = subprocess.run(cmd, text=True, capture_output=True)
        if completed.returncode == 0:
            raise AssertionError(
                "BiWFA accepted negative-match ends-free configuration:\n"
                + " ".join(cmd)
            )
        if completed.stdout != "" or completed.stderr != NEGATIVE_MATCH_MESSAGE:
            raise AssertionError(
                "BiWFA rejected negative-match ends-free with an unexpected message:\n"
                + " ".join(cmd)
                + f"\n\nstdout: {completed.stdout!r}\nstderr: {completed.stderr!r}"
            )


def check_heuristic_smoke(align_benchmark: Path, workdir: Path) -> None:
    free_ends = (2, 2, 2, 2)
    span = "ends-free,%d,%d,%d,%d" % free_ends
    cases = [
        Case("TTACGTACGTG", "AACGTTCGTGGG"),
        Case("CCGATTACA", "TGATTTACAAA"),
    ]
    metrics = [
        metric for metric in METRICS
        if metric[0] in ("affine", "affine2p")
    ]
    for heuristic_label, heuristic_args in HEURISTIC_SMOKE_CONFIGS:
        check_cases(
            align_benchmark, workdir, f"heuristic.{heuristic_label}",
            free_ends, span, cases, metrics=metrics, extra_args=heuristic_args)


def check_recursive_smoke(align_benchmark: Path, workdir: Path) -> None:
    pattern, text = divergent_pair(900)
    cases = [Case("TTTT" + pattern + "GGGG", "CCCC" + text + "AAAA")]
    input_path = workdir / "recursive.seq"
    output_path = workdir / "recursive.alg"
    write_cases(input_path, cases)
    cmd = [
        str(align_benchmark),
        "--input",
        str(input_path),
        "--output-full",
        str(output_path),
        "--algorithm",
        "gap-affine-wfa",
        "--wfa-memory=ultralow",
        "--wfa-span=ends-free,4,4,4,4",
        "--verbose=3",
    ]
    completed = subprocess.run(cmd, text=True, capture_output=True)
    if completed.returncode != 0:
        raise RuntimeError("recursive smoke failed:\n" + " ".join(cmd))
    debug = completed.stdout + completed.stderr
    if "[WFA::BiAlign][Recursion=1]" not in debug:
        raise AssertionError(
            "long BiWFA case did not show recursive divide-and-conquer\n"
            + " ".join(cmd)
        )


def main(argv: list[str]) -> int:
    align_benchmark = find_align_benchmark(argv)
    with tempfile.TemporaryDirectory() as tmpdir:
        workdir = Path(tmpdir)
        label, free_ends = GLOBAL_CONFIG
        check_config(align_benchmark, workdir, label, free_ends, "global")
        for label, free_ends in FREE_END_CONFIGS:
            span = "ends-free,%d,%d,%d,%d" % free_ends
            check_config(align_benchmark, workdir, label, free_ends, span)
        check_tiny_exhaustive(align_benchmark, workdir)
        check_negative_match_rejection(align_benchmark, workdir)
        check_heuristic_smoke(align_benchmark, workdir)
        check_recursive_smoke(align_benchmark, workdir)
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv))
