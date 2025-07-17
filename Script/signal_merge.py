#!/usr/bin/env python3
from pathlib import Path
import argparse, sys
import pod5 as p5
from itertools import zip_longest
from tqdm import tqdm

def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description="Interleave-merge POD5 files: nth read from each input → output")
    p.add_argument("inputs", metavar="INPUT.pod5", nargs="+",
                   help="Two or more input POD5 files produced in the SAME order")
    p.add_argument("-o", "--output", required=True, metavar="OUTPUT.pod5",
                   help="Destination POD5 file (overwritten if it exists)")
    p.add_argument("--progress", action="store_true",
                   help="Show a tqdm progress bar")
    return p.parse_args()

def open_readers(paths):
    readers = []
    for path in paths:
        if not Path(path).is_file():
            sys.exit(f"[ERROR] {path} does not exist or is not a file")
        reader = p5.Reader(path)
        readers.append(reader)
    return readers

def main():
    args = parse_args()
    if len(args.inputs) < 2:
        sys.exit("[ERROR] Provide at least TWO input POD5 files")

    readers = open_readers(args.inputs)
    gens     = [r.reads() for r in readers]          # generators of ReadRecord ✔
    try:
        total_reads = readers[0].num_reads           # v≥0.3.20
    except AttributeError:
        total_reads = None

    bar = tqdm(total=total_reads, disable=not args.progress, unit="variant")

    with p5.Writer(args.output) as writer:
        for idx, read_tuple in enumerate(zip_longest(*gens, fillvalue=None)):
            if any(r is None for r in read_tuple):
                sys.exit(f"[ERROR] Input files differ in length (stopped at read {idx}).")

            # --- FIXED LINE ---------------------------------------------------
            writer.add_reads([rec.to_read() for rec in read_tuple])
            # ------------------------------------------------------------------

            if args.progress:
                bar.update()

    for r in readers:
        r.close()
    if args.progress:
        bar.close()

    print(f"[✓] Wrote interleaved dataset to: {args.output}")

if __name__ == "__main__":
    main()