#!/usr/bin/env python3
"""Singletrack correctness checks."""

from __future__ import annotations

from dataclasses import dataclass
import random
import re
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path


DNA = "ACGT"
TINY_DNA = "AC"
NEG_INF = -10**12
CIGAR_RE = re.compile(r"(\d*)([MXID])")

FREE_END_CONFIGS = [
    ("all", (4, 4, 4, 4)),
    ("pattern-begin", (4, 0, 0, 0)),
    ("pattern-end", (0, 4, 0, 0)),
    ("text-begin", (0, 0, 4, 0)),
    ("text-end", (0, 0, 0, 4)),
    ("pattern-only", (4, 4, 0, 0)),
    ("text-only", (0, 0, 4, 4)),
]


PENALTY_CONFIGS = [
    (
        "affine.default",
        "gap-affine-wfa",
        [],
        {"kind": "affine", "match": 0, "mismatch": 4, "open": 6, "extend": 2},
    ),
    (
        "affine.negmatch",
        "gap-affine-wfa",
        ["--affine-penalties=-5,1,2,1"],
        {"kind": "affine", "match": -5, "mismatch": 1, "open": 2, "extend": 1},
    ),
    (
        "affine2p.default",
        "gap-affine2p-wfa",
        [],
        {
            "kind": "affine2p",
            "match": 0,
            "mismatch": 4,
            "open1": 6,
            "extend1": 2,
            "open2": 24,
            "extend2": 1,
        },
    ),
    (
        "affine2p.negmatch",
        "gap-affine2p-wfa",
        ["--affine2p-penalties=-5,1,2,1,5,1"],
        {
            "kind": "affine2p",
            "match": -5,
            "mismatch": 1,
            "open1": 2,
            "extend1": 1,
            "open2": 5,
            "extend2": 1,
        },
    ),
]


NO_HEURISTIC_CONFIG = ("none", [])

SUPPORTED_HEURISTIC_CONFIGS = [
    (
        "wfadaptive",
        ["--wfa-heuristic=wfa-adaptive", "--wfa-heuristic-parameters=10,50,1"],
    ),
    (
        "wfmash",
        ["--wfa-heuristic=wfmash", "--wfa-heuristic-parameters=10,50,1"],
    ),
    (
        "xdrop",
        ["--wfa-heuristic=xdrop", "--wfa-heuristic-parameters=100,1"],
    ),
    (
        "zdrop",
        ["--wfa-heuristic=zdrop", "--wfa-heuristic-parameters=100,1"],
    ),
]

UNSUPPORTED_SINGLETRACK_HEURISTIC_CONFIGS = [
    (
        "banded-static",
        ["--wfa-heuristic=banded-static", "--wfa-heuristic-parameters=-10,10"],
    ),
    (
        "banded-adaptive",
        ["--wfa-heuristic=banded-adaptive", "--wfa-heuristic-parameters=-10,10,1"],
    ),
]


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


def find_align_benchmark(argv: list[str]) -> Path:
    if len(argv) > 1:
        candidate = Path(argv[1])
        if candidate.is_dir():
            candidate = candidate / "align_benchmark"
        if candidate.is_file():
            return candidate
    script_dir = Path(__file__).resolve().parent
    root = script_dir.parent
    for candidate in (root / "bin" / "align_benchmark", root / "build" / "align_benchmark"):
        if candidate.is_file():
            return candidate
    raise SystemExit("[Error] Binaries not built. Please run cmake or make")


def random_dna(rng: random.Random, length: int) -> str:
    return "".join(rng.choice(DNA) for _ in range(length))


def enumerate_sequences(alphabet: str, max_length: int) -> list[str]:
    sequences: list[str] = []
    current = [""]
    for _ in range(max_length):
        current = [prefix + base for prefix in current for base in alphabet]
        sequences.extend(current)
    return sequences


def mutate(rng: random.Random, sequence: str, max_edits: int) -> str:
    bases = list(sequence)
    for _ in range(rng.randint(0, max_edits)):
        op = rng.choice("XXIIDM") if bases else "I"
        if op == "X":
            idx = rng.randrange(len(bases))
            bases[idx] = rng.choice([base for base in DNA if base != bases[idx]])
        elif op == "I":
            bases.insert(rng.randrange(len(bases) + 1), rng.choice(DNA))
        elif op == "D":
            del bases[rng.randrange(len(bases))]
    return "".join(bases) or rng.choice(DNA)


def exhaustive_tiny_pairs(free_ends: tuple[int, int, int, int] | None = None) -> list[tuple[str, str]]:
    sequences = enumerate_sequences(TINY_DNA, 5)
    pairs: list[tuple[str, str]] = []
    for pattern in sequences:
        for text in sequences:
            if free_ends is not None:
                p0, pf, t0, tf = free_ends
                if p0 > len(pattern) or pf > len(pattern):
                    continue
                if t0 > len(text) or tf > len(text):
                    continue
            pairs.append((pattern, text))
    return pairs


def endsfree_adversarial_pairs(free_ends: tuple[int, int, int, int]) -> list[tuple[str, str]]:
    p0, pf, t0, tf = free_ends
    core = "ACGTACGTAC"
    pairs = [
        ("AAAAAAA", "AAAATAAA"),
        ("ACACACAC", "ACACGTAC"),
        ("TTTTAAAACCCC", "AAAACCCCGGGG"),
        ("GGGGACGTACGT", "ACGTACGTCCCC"),
    ]
    if p0 > 0:
        pairs.append(("T" * p0 + core, core))
        pairs.append(("T" * (p0 + 1) + core, core))
    if pf > 0:
        pairs.append((core + "G" * pf, core))
        pairs.append((core + "G" * (pf + 1), core))
    if t0 > 0:
        pairs.append((core, "T" * t0 + core))
        pairs.append((core, "T" * (t0 + 1) + core))
    if tf > 0:
        pairs.append((core, core + "G" * tf))
        pairs.append((core, core + "G" * (tf + 1)))
    if p0 > 0 and tf > 0:
        pairs.append(("C" * p0 + core, core + "G" * tf))
        pairs.append(("C" * (p0 + 1) + core, core + "G" * (tf + 1)))
    if t0 > 0 and pf > 0:
        pairs.append((core + "C" * pf, "G" * t0 + core))
        pairs.append((core + "C" * (pf + 1), "G" * (t0 + 1) + core))
    return pairs


def endsfree_pairs(
    free_ends: tuple[int, int, int, int],
    include_exhaustive: bool = True,
) -> list[tuple[str, str]]:
    rng = random.Random(0x5A17)
    pairs = [
        ("AAAACCCCGGGG", "TTTTAAAACCCCGGGGAAAA"),
        ("AAAACCCCGGGGTTTT", "GGGGAAAACCCCGGGG"),
        (
            "CGTATTTACTTGGGACGCATCGTATGCACACATGCATAGAGTTTTCCCACGCGCCGGACCTGAATGACA",
            "CAGCTTACTTGCGACGCTCGTGTGCACAATGCATAGAGTTTTCCACGCGCCGAGACCTGTAACG",
        ),
    ]
    pairs.extend(endsfree_adversarial_pairs(free_ends))
    for _ in range(40):
        base = random_dna(rng, rng.randint(4, 28))
        mate = mutate(rng, base, 6)
        pairs.append(
            (
                random_dna(rng, rng.randint(0, free_ends[0]))
                + base
                + random_dna(rng, rng.randint(0, free_ends[1])),
                random_dna(rng, rng.randint(0, free_ends[2]))
                + mate
                + random_dna(rng, rng.randint(0, free_ends[3])),
            )
        )
    if include_exhaustive:
        pairs.extend(exhaustive_tiny_pairs(free_ends))
    return pairs


def extension_adversarial_pairs() -> list[tuple[str, str]]:
    return [
        ("AAAAAAAATTTT", "AAAATAAAGGGG"),
        ("ACACACACGGGG", "ACACGTACCCCC"),
        ("ACGTACGT" + "T" * 8, "ACGTACGT" + "G" * 8),
        ("ACGTACGTGGGG", "ACGTTACGTCCCC"),
        ("TTTTAAAACCCCGGGG", "TTTTAAAAGGGGCCCC"),
        ("GATTACAGGGGGG", "GATTTACACCCCCC"),
    ]


def extension_pairs(include_exhaustive: bool = True) -> list[tuple[str, str]]:
    rng = random.Random(0xE617)
    pairs = [
        ("AAAACCCCGGGGTTTT", "AAAACCCCGGGGAAAA"),
        ("ACGTACGTGGGG", "ACGTACGTTTTT"),
        (
            "GGAGTCAGGACGAACGATATGTCGCTGCACCTATTCACTTGAGTTGCCAGTATAAATGAACGATGAGCTACCGGTGCTACACAAATCCAATTGCGAG",
            "GGAGTCAGACGAACCGGATATGTCGTGCAACCTATTCACTATGAGTTGACCGTATAAATGAACGATAGAGCTACGGTGCTGTATTAACATACTCGTGAGGCAAAG",
        ),
    ]
    pairs.extend(extension_adversarial_pairs())
    for _ in range(40):
        base = random_dna(rng, rng.randint(4, 32))
        mate = mutate(rng, base, 6)
        pairs.append(
            (
                base + random_dna(rng, rng.randint(0, 10)),
                mate + random_dna(rng, rng.randint(0, 10)),
            )
        )
    if include_exhaustive:
        pairs.extend(exhaustive_tiny_pairs())
    return pairs


def write_pairs(path: Path, pairs: list[tuple[str, str]]) -> None:
    with path.open("w") as handle:
        for pattern, text in pairs:
            print(f">{pattern}", file=handle)
            print(f"<{text}", file=handle)


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


def validate_cigar(
    pattern: str,
    text: str,
    cigar: str,
    allow_partial: bool,
) -> list[tuple[int, int]]:
    v = 0
    h = 0
    coords = [(0, 0)]
    for op in parse_cigar(cigar):
        if op == "M":
            if v >= len(pattern) or h >= len(text) or pattern[v] != text[h]:
                raise AssertionError(f"invalid M at ({v},{h}) in {cigar}")
            v += 1
            h += 1
        elif op == "X":
            if v >= len(pattern) or h >= len(text) or pattern[v] == text[h]:
                raise AssertionError(f"invalid X at ({v},{h}) in {cigar}")
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
    if not allow_partial and (v != len(pattern) or h != len(text)):
        raise AssertionError(
            f"CIGAR ends at ({v},{h}), expected ({len(pattern)},{len(text)}) in {cigar}"
        )
    return coords


def match_score(penalties: dict[str, int], extension: bool) -> int:
    if extension and penalties["match"] == 0:
        return 1
    return -penalties["match"]


def score_segment(
    ops: list[str],
    penalties: dict[str, int],
    extension: bool = False,
) -> int:
    score = 0
    if penalties["kind"] == "affine":
        last = ""
        for op in ops:
            if op == "M":
                score += match_score(penalties, extension)
            elif op == "X":
                score -= penalties["mismatch"]
            elif op in ("I", "D"):
                score -= penalties["extend"]
                if last != op:
                    score -= penalties["open"]
            last = op
        return score
    idx = 0
    while idx < len(ops):
        op = ops[idx]
        end = idx + 1
        while end < len(ops) and ops[end] == op:
            end += 1
        length = end - idx
        if op == "M":
            score += match_score(penalties, extension) * length
        elif op == "X":
            score -= penalties["mismatch"] * length
        elif op in ("I", "D"):
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
    coords = validate_cigar(pattern, text, cigar, allow_partial=False)
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
        score_segment(ops[start:end], penalties)
        for start in starts
        for end in ends
        if end >= start
    ]
    if not candidates:
        raise AssertionError(f"CIGAR has no legal ends-free segment: {cigar}")
    return max(candidates)


def cigar_extension_score(pattern: str, text: str, cigar: str, penalties: dict[str, int]) -> int:
    if cigar == "-":
        return 0
    validate_cigar(pattern, text, cigar, allow_partial=True)
    return score_segment(parse_cigar(cigar), penalties, extension=True)


def dp_endsfree(
    pattern: str,
    text: str,
    penalties: dict[str, int],
    free_ends: tuple[int, int, int, int],
) -> int:
    plen = len(pattern)
    tlen = len(text)
    p0, pf, t0, tf = free_ends
    if penalties["kind"] == "affine":
        m = [[NEG_INF] * (plen + 1) for _ in range(tlen + 1)]
        d = [[NEG_INF] * (plen + 1) for _ in range(tlen + 1)]
        i = [[NEG_INF] * (plen + 1) for _ in range(tlen + 1)]
        m[0][0] = 0
        for v in range(1, plen + 1):
            if v <= p0:
                m[0][v] = 0
            else:
                d[0][v] = max(
                    m[0][v - 1] - penalties["open"] - penalties["extend"],
                    d[0][v - 1] - penalties["extend"],
                )
                m[0][v] = d[0][v]
        for h in range(1, tlen + 1):
            if h <= t0:
                m[h][0] = 0
            else:
                i[h][0] = max(
                    m[h - 1][0] - penalties["open"] - penalties["extend"],
                    i[h - 1][0] - penalties["extend"],
                )
                m[h][0] = i[h][0]
        for h in range(1, tlen + 1):
            for v in range(1, plen + 1):
                d[h][v] = max(
                    m[h][v - 1] - penalties["open"] - penalties["extend"],
                    d[h][v - 1] - penalties["extend"],
                )
                i[h][v] = max(
                    m[h - 1][v] - penalties["open"] - penalties["extend"],
                    i[h - 1][v] - penalties["extend"],
                )
                diag = m[h - 1][v - 1]
                diag += match_score(penalties, extension=False) if pattern[v - 1] == text[h - 1] else -penalties["mismatch"]
                m[h][v] = max(diag, d[h][v], i[h][v])
        return max(
            [m[h][plen] for h in range(max(0, tlen - tf), tlen + 1)]
            + [m[tlen][v] for v in range(max(0, plen - pf), plen + 1)]
        )

    states = {
        name: [[NEG_INF] * (plen + 1) for _ in range(tlen + 1)]
        for name in ("M", "D1", "D2", "I1", "I2")
    }
    m = states["M"]
    m[0][0] = 0
    for v in range(1, plen + 1):
        if v <= p0:
            m[0][v] = 0
        else:
            states["D1"][0][v] = max(
                m[0][v - 1] - penalties["open1"] - penalties["extend1"],
                states["D1"][0][v - 1] - penalties["extend1"],
            )
            states["D2"][0][v] = max(
                m[0][v - 1] - penalties["open2"] - penalties["extend2"],
                states["D2"][0][v - 1] - penalties["extend2"],
            )
            m[0][v] = max(states["D1"][0][v], states["D2"][0][v])
    for h in range(1, tlen + 1):
        if h <= t0:
            m[h][0] = 0
        else:
            states["I1"][h][0] = max(
                m[h - 1][0] - penalties["open1"] - penalties["extend1"],
                states["I1"][h - 1][0] - penalties["extend1"],
            )
            states["I2"][h][0] = max(
                m[h - 1][0] - penalties["open2"] - penalties["extend2"],
                states["I2"][h - 1][0] - penalties["extend2"],
            )
            m[h][0] = max(states["I1"][h][0], states["I2"][h][0])
    for h in range(1, tlen + 1):
        for v in range(1, plen + 1):
            states["D1"][h][v] = max(
                m[h][v - 1] - penalties["open1"] - penalties["extend1"],
                states["D1"][h][v - 1] - penalties["extend1"],
            )
            states["D2"][h][v] = max(
                m[h][v - 1] - penalties["open2"] - penalties["extend2"],
                states["D2"][h][v - 1] - penalties["extend2"],
            )
            states["I1"][h][v] = max(
                m[h - 1][v] - penalties["open1"] - penalties["extend1"],
                states["I1"][h - 1][v] - penalties["extend1"],
            )
            states["I2"][h][v] = max(
                m[h - 1][v] - penalties["open2"] - penalties["extend2"],
                states["I2"][h - 1][v] - penalties["extend2"],
            )
            diag = m[h - 1][v - 1]
            diag += match_score(penalties, extension=False) if pattern[v - 1] == text[h - 1] else -penalties["mismatch"]
            m[h][v] = max(
                diag,
                states["D1"][h][v],
                states["D2"][h][v],
                states["I1"][h][v],
                states["I2"][h][v],
            )
    return max(
        [m[h][plen] for h in range(max(0, tlen - tf), tlen + 1)]
        + [m[tlen][v] for v in range(max(0, plen - pf), plen + 1)]
    )


def dp_extension(pattern: str, text: str, penalties: dict[str, int]) -> int:
    plen = len(pattern)
    tlen = len(text)
    if penalties["kind"] == "affine":
        m = [[NEG_INF] * (plen + 1) for _ in range(tlen + 1)]
        d = [[NEG_INF] * (plen + 1) for _ in range(tlen + 1)]
        i = [[NEG_INF] * (plen + 1) for _ in range(tlen + 1)]
        m[0][0] = 0
        best = 0
        for v in range(1, plen + 1):
            d[0][v] = max(
                m[0][v - 1] - penalties["open"] - penalties["extend"],
                d[0][v - 1] - penalties["extend"],
            )
            m[0][v] = d[0][v]
            best = max(best, m[0][v])
        for h in range(1, tlen + 1):
            i[h][0] = max(
                m[h - 1][0] - penalties["open"] - penalties["extend"],
                i[h - 1][0] - penalties["extend"],
            )
            m[h][0] = i[h][0]
            best = max(best, m[h][0])
        for h in range(1, tlen + 1):
            for v in range(1, plen + 1):
                d[h][v] = max(
                    m[h][v - 1] - penalties["open"] - penalties["extend"],
                    d[h][v - 1] - penalties["extend"],
                )
                i[h][v] = max(
                    m[h - 1][v] - penalties["open"] - penalties["extend"],
                    i[h - 1][v] - penalties["extend"],
                )
                diag = m[h - 1][v - 1]
                diag += match_score(penalties, extension=True) if pattern[v - 1] == text[h - 1] else -penalties["mismatch"]
                m[h][v] = max(diag, d[h][v], i[h][v])
                best = max(best, m[h][v])
        return best

    states = {
        name: [[NEG_INF] * (plen + 1) for _ in range(tlen + 1)]
        for name in ("M", "D1", "D2", "I1", "I2")
    }
    m = states["M"]
    m[0][0] = 0
    best = 0
    for v in range(1, plen + 1):
        states["D1"][0][v] = max(
            m[0][v - 1] - penalties["open1"] - penalties["extend1"],
            states["D1"][0][v - 1] - penalties["extend1"],
        )
        states["D2"][0][v] = max(
            m[0][v - 1] - penalties["open2"] - penalties["extend2"],
            states["D2"][0][v - 1] - penalties["extend2"],
        )
        m[0][v] = max(states["D1"][0][v], states["D2"][0][v])
        best = max(best, m[0][v])
    for h in range(1, tlen + 1):
        states["I1"][h][0] = max(
            m[h - 1][0] - penalties["open1"] - penalties["extend1"],
            states["I1"][h - 1][0] - penalties["extend1"],
        )
        states["I2"][h][0] = max(
            m[h - 1][0] - penalties["open2"] - penalties["extend2"],
            states["I2"][h - 1][0] - penalties["extend2"],
        )
        m[h][0] = max(states["I1"][h][0], states["I2"][h][0])
        best = max(best, m[h][0])
    for h in range(1, tlen + 1):
        for v in range(1, plen + 1):
            states["D1"][h][v] = max(
                m[h][v - 1] - penalties["open1"] - penalties["extend1"],
                states["D1"][h][v - 1] - penalties["extend1"],
            )
            states["D2"][h][v] = max(
                m[h][v - 1] - penalties["open2"] - penalties["extend2"],
                states["D2"][h][v - 1] - penalties["extend2"],
            )
            states["I1"][h][v] = max(
                m[h - 1][v] - penalties["open1"] - penalties["extend1"],
                states["I1"][h - 1][v] - penalties["extend1"],
            )
            states["I2"][h][v] = max(
                m[h - 1][v] - penalties["open2"] - penalties["extend2"],
                states["I2"][h - 1][v] - penalties["extend2"],
            )
            diag = m[h - 1][v - 1]
            diag += match_score(penalties, extension=True) if pattern[v - 1] == text[h - 1] else -penalties["mismatch"]
            m[h][v] = max(
                diag,
                states["D1"][h][v],
                states["D2"][h][v],
                states["I1"][h][v],
                states["I2"][h][v],
            )
            best = max(best, m[h][v])
    return best


def run_align_benchmark(
    align_benchmark: Path,
    workdir: Path,
    name: str,
    pairs: list[tuple[str, str]],
    algorithm: str,
    memory_mode: str,
    span: str,
    extra_args: list[str],
) -> RunResult:
    input_path = workdir / f"{name}.seq"
    output_path = workdir / f"{name}.{memory_mode}.alg"
    write_pairs(input_path, pairs)
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
        *extra_args,
    ]
    completed = subprocess.run(cmd, cwd=align_benchmark.parents[1], text=True, capture_output=True)
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
    if len(rows) != len(pairs):
        raise AssertionError(f"expected {len(pairs)} outputs, got {len(rows)}")
    return RunResult(
        rows=rows,
        command=cmd,
        input_path=input_path,
        output_path=output_path,
    )


def heuristic_configs_for_penalties(
    penalties: dict[str, int],
) -> list[tuple[str, list[str]]]:
    if penalties["match"] == 0:
        return [NO_HEURISTIC_CONFIG, *SUPPORTED_HEURISTIC_CONFIGS]
    return [NO_HEURISTIC_CONFIG]


def reject_singletrack_config(
    align_benchmark: Path,
    workdir: Path,
    name: str,
    extra_args: list[str],
) -> None:
    pairs = [("ACGTACGT", "ACGTTACGT")]
    input_path = workdir / f"reject.{name}.seq"
    output_path = workdir / f"reject.{name}.alg"
    write_pairs(input_path, pairs)
    cmd = [
        str(align_benchmark),
        "--input",
        str(input_path),
        "--output-full",
        str(output_path),
        "--algorithm",
        "gap-affine-wfa",
        "--wfa-memory=singletrack",
        *extra_args,
    ]
    completed = subprocess.run(cmd, cwd=align_benchmark.parents[1], text=True, capture_output=True)
    if completed.returncode == 0:
        raise AssertionError(
            "Singletrack accepted unsupported heuristic configuration:\n"
            + " ".join(cmd)
        )
    message = completed.stdout + completed.stderr
    if "Singletrack" not in message or "banded" not in message:
        raise AssertionError(
            "Singletrack rejected unsupported heuristic with an unexpected message:\n"
            + " ".join(cmd)
            + "\n\n"
            + message
        )


def fail_with_replay(
    message: str,
    high: RunResult,
    singletrack: RunResult,
    case_index: int,
    pattern: str,
    text: str,
    high_cigar: str,
    singletrack_cigar: str,
    observed: int,
    expected: int,
) -> None:
    raise AssertionError(
        f"{message}\n"
        f"  case: {case_index}\n"
        f"  observed: {observed}\n"
        f"  expected: {expected}\n"
        f"  pattern: {pattern}\n"
        f"  text: {text}\n"
        f"  high_cigar: {high_cigar}\n"
        f"  singletrack_cigar: {singletrack_cigar}\n"
        f"  high_input: {high.input_path}\n"
        f"  high_output: {high.output_path}\n"
        f"  singletrack_input: {singletrack.input_path}\n"
        f"  singletrack_output: {singletrack.output_path}\n"
        f"  replay_high: {' '.join(high.command)}\n"
        f"  replay_singletrack: {' '.join(singletrack.command)}"
    )


def skip_unreachable_row(
    message: str,
    high: RunResult,
    singletrack: RunResult,
    case_index: int,
    high_row: AlignmentRow,
    singletrack_row: AlignmentRow,
) -> bool:
    if (
        high_row.pattern != singletrack_row.pattern
        or high_row.text != singletrack_row.text
    ):
        fail_with_replay(
            f"{message} case {case_index} output sequence differs",
            high,
            singletrack,
            case_index,
            high_row.pattern,
            high_row.text,
            high_row.cigar,
            singletrack_row.cigar,
            singletrack_row.score,
            high_row.score,
        )
    if high_row.cigar == "-" or singletrack_row.cigar == "-":
        if (
            high_row.score != -1
            or singletrack_row.score != -1
            or high_row.cigar != singletrack_row.cigar
        ):
            fail_with_replay(
                f"{message} case {case_index} unreachable output differs",
                high,
                singletrack,
                case_index,
                high_row.pattern,
                high_row.text,
                high_row.cigar,
                singletrack_row.cigar,
                singletrack_row.score,
                high_row.score,
            )
        return True
    return False


def check_endsfree(align_benchmark: Path, workdir: Path) -> None:
    for free_label, free_ends in FREE_END_CONFIGS:
        span = "ends-free,%d,%d,%d,%d" % free_ends
        for label, algorithm, extra_args, penalties in PENALTY_CONFIGS:
            for heuristic_label, heuristic_args in heuristic_configs_for_penalties(penalties):
                include_exhaustive = heuristic_label == "none"
                pairs = endsfree_pairs(free_ends, include_exhaustive=include_exhaustive)
                config_label = (
                    label if heuristic_label == "none" else f"{label}.{heuristic_label}"
                )
                combined_args = [*extra_args, *heuristic_args]
                print(f">>> Proving Singletrack ends-free {config_label} {free_label}")
                high = run_align_benchmark(
                    align_benchmark,
                    workdir,
                    f"endsfree.{free_label}.{config_label}",
                    pairs,
                    algorithm,
                    "high",
                    span,
                    combined_args,
                )
                singletrack = run_align_benchmark(
                    align_benchmark,
                    workdir,
                    f"endsfree.{free_label}.{config_label}",
                    pairs,
                    algorithm,
                    "singletrack",
                    span,
                    combined_args,
                )
                for idx, (high_row, singletrack_row) in enumerate(
                    zip(high.rows, singletrack.rows), start=1
                ):
                    if skip_unreachable_row(
                        f"ends-free {config_label} {free_label}",
                        high,
                        singletrack,
                        idx,
                        high_row,
                        singletrack_row,
                    ):
                        continue
                    pattern = high_row.pattern
                    text = high_row.text
                    high_cigar = high_row.cigar
                    singletrack_cigar = singletrack_row.cigar
                    baseline = cigar_endsfree_score(pattern, text, high_cigar, penalties, free_ends)
                    observed = cigar_endsfree_score(
                        pattern, text, singletrack_cigar, penalties, free_ends
                    )
                    if observed != baseline:
                        fail_with_replay(
                            f"ends-free {config_label} {free_label} case {idx} differs from "
                            "high-memory WFA",
                            high,
                            singletrack,
                            idx,
                            pattern,
                            text,
                            high_cigar,
                            singletrack_cigar,
                            observed,
                            baseline,
                        )
                    if heuristic_label == "none" and penalties["match"] == 0:
                        expected = dp_endsfree(pattern, text, penalties, free_ends)
                        if observed != expected:
                            fail_with_replay(
                                f"ends-free {config_label} {free_label} case {idx} "
                                "failed DP oracle",
                                high,
                                singletrack,
                                idx,
                                pattern,
                                text,
                                high_cigar,
                                singletrack_cigar,
                                observed,
                                expected,
                            )


def check_extension(align_benchmark: Path, workdir: Path) -> None:
    for label, algorithm, extra_args, penalties in PENALTY_CONFIGS:
        for heuristic_label, heuristic_args in heuristic_configs_for_penalties(penalties):
            include_exhaustive = penalties["match"] == 0 and heuristic_label == "none"
            pairs = extension_pairs(include_exhaustive=include_exhaustive)
            config_label = label if heuristic_label == "none" else f"{label}.{heuristic_label}"
            combined_args = [*extra_args, *heuristic_args]
            print(f">>> Proving Singletrack extension {config_label}")
            high = run_align_benchmark(
                align_benchmark,
                workdir,
                f"extension.{config_label}",
                pairs,
                algorithm,
                "high",
                "extension",
                combined_args,
            )
            singletrack = run_align_benchmark(
                align_benchmark,
                workdir,
                f"extension.{config_label}",
                pairs,
                algorithm,
                "singletrack",
                "extension",
                combined_args,
            )
            for idx, (high_row, singletrack_row) in enumerate(
                zip(high.rows, singletrack.rows), start=1
            ):
                if skip_unreachable_row(
                    f"extension {config_label}",
                    high,
                    singletrack,
                    idx,
                    high_row,
                    singletrack_row,
                ):
                    continue
                pattern = high_row.pattern
                text = high_row.text
                high_cigar = high_row.cigar
                singletrack_cigar = singletrack_row.cigar
                baseline = cigar_extension_score(pattern, text, high_cigar, penalties)
                observed = cigar_extension_score(pattern, text, singletrack_cigar, penalties)
                if observed != baseline:
                    fail_with_replay(
                        f"extension {config_label} case {idx} differs from high-memory WFA",
                        high,
                        singletrack,
                        idx,
                        pattern,
                        text,
                        high_cigar,
                        singletrack_cigar,
                        observed,
                        baseline,
                    )
                if heuristic_label == "none" and penalties["match"] == 0:
                    expected = dp_extension(pattern, text, penalties)
                    if observed != expected:
                        fail_with_replay(
                            f"extension {config_label} case {idx} failed DP oracle",
                            high,
                            singletrack,
                            idx,
                            pattern,
                            text,
                            high_cigar,
                            singletrack_cigar,
                            observed,
                            expected,
                        )


def check_unsupported_heuristics(align_benchmark: Path, workdir: Path) -> None:
    for label, heuristic_args in UNSUPPORTED_SINGLETRACK_HEURISTIC_CONFIGS:
        print(f">>> Proving Singletrack rejects {label}")
        reject_singletrack_config(align_benchmark, workdir, label, heuristic_args)


def main(argv: list[str]) -> int:
    align_benchmark = find_align_benchmark(argv).resolve()
    workdir = Path(tempfile.mkdtemp(prefix="wfa.singletrack.proof."))
    try:
        check_endsfree(align_benchmark, workdir)
        check_extension(align_benchmark, workdir)
        check_unsupported_heuristics(align_benchmark, workdir)
        print(f">>> Singletrack proof tests passed ({workdir})")
    except Exception:
        print(f">>> Singletrack proof temp files preserved at {workdir}", file=sys.stderr)
        raise
    else:
        shutil.rmtree(workdir)
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv))
