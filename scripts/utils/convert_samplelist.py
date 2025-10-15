#!/usr/bin/env python3
import argparse
import os
import re
import sys
from pathlib import Path

PREFER_DIRS = [
    'AP', 'scRNA-seq', 'WES', 'ATAC', 'Raw_ATAC', 'ATAC2'
]

FASTQ_RE = re.compile(r"\.(?:fastq|fq)(?:\.gz)?$", re.IGNORECASE)


def parse_args():
    p = argparse.ArgumentParser(description='Convert legacy samplelist to 3-column samples.list (sample_id\tR1\tR2).')
    p.add_argument('--in', dest='input', required=True, help='Path to legacy samplelist file.')
    p.add_argument('--out', dest='output', required=True, help='Path to write standardized samples.list.')
    p.add_argument('--search-root', dest='search_root', default='.', help='Root directory to search for FASTQ files.')
    return p.parse_args()


def prefer_score(path: str) -> int:
    # Lower score is better; prefer paths containing earlier entries in PREFER_DIRS
    score = 1000
    for i, d in enumerate(PREFER_DIRS):
        if f"/{d}/" in path or path.startswith(d + "/"):
            score = min(score, i)
    return score


def find_file(search_root: Path, token: str) -> str | None:
    p = Path(token)
    if p.is_file():
        return str(p)
    rp = (search_root / token)
    if rp.is_file():
        return str(rp)
    basename = os.path.basename(token)
    candidates = []
    for root, _, files in os.walk(search_root):
        if basename in files:
            candidates.append(os.path.join(root, basename))
    if not candidates:
        return None
    # Sort by prefer score then by path length
    candidates.sort(key=lambda x: (prefer_score(x), len(x)))
    return candidates[0]


def extract_fastqs(tokens: list[str]) -> tuple[str | None, str | None]:
    fastqs = [t for t in tokens if FASTQ_RE.search(t)]
    r1 = next((t for t in fastqs if 'R1' in t), None)
    r2 = next((t for t in fastqs if 'R2' in t), None)
    return r1, r2


def convert(input_path: Path, output_path: Path, search_root: Path) -> int:
    lines_out = []
    missing = 0
    with input_path.open('r', encoding='utf-8', errors='ignore') as f:
        for line in f:
            raw = line.strip()
            if not raw or raw.startswith('#'):
                continue
            tokens = re.split(r"\s+", raw)
            if not tokens:
                continue
            sample_id = tokens[0]
            r1_token, r2_token = extract_fastqs(tokens)
            r1_path = find_file(search_root, r1_token) if r1_token else None
            r2_path = find_file(search_root, r2_token) if r2_token else None

            if r1_path:
                r1_path = os.path.relpath(r1_path, search_root)
            if r2_path:
                r2_path = os.path.relpath(r2_path, search_root)

            if not (r1_path and r2_path):
                missing += 1
                print(f"[WARN] Missing R1/R2 for sample {sample_id}: R1={r1_token} -> {r1_path}, R2={r2_token} -> {r2_path}", file=sys.stderr)
            lines_out.append((sample_id, r1_path or (r1_token or ''), r2_path or (r2_token or '')))

    output_path.parent.mkdir(parents=True, exist_ok=True)
    with output_path.open('w', encoding='utf-8') as out:
        out.write('# sample_id\tR1_fastq_path\tR2_fastq_path\n')
        for sid, r1, r2 in lines_out:
            out.write(f"{sid}\t{r1}\t{r2}\n")

    return missing


def main():
    args = parse_args()
    input_path = Path(args.input).resolve()
    output_path = Path(args.output)
    search_root = Path(args.search_root).resolve()

    if not input_path.exists():
        print(f"[ERROR] Input samplelist not found: {input_path}", file=sys.stderr)
        sys.exit(2)
    missing = convert(input_path, output_path, search_root)
    print(f"[INFO] Wrote standardized samples.list to: {output_path}")
    if missing:
        print(f"[INFO] {missing} samples had unresolved paths; please adjust manually.")

if __name__ == '__main__':
    main()
