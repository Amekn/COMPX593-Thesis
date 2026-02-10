#!/usr/bin/env python3
"""
pod5Merger.py

combine multiple POD5 files into a single POD5 file.

Requires:
  pip install pod5

Example:
  python pod5Merger.py <output_pod5> <input_pod5_1> <input_pod5_2> ... <input_pod5_n>
"""

from __future__ import annotations

import argparse
from pathlib import Path
from typing import List, Optional

import pod5 as p5

def build_parser() -> argparse.ArgumentParser:
    ap = argparse.ArgumentParser(
        description="combine multiple POD5 files into a single POD5 file.")
    ap.add_argument("output_pod5", type=Path, help="Output .pod5 file")
    ap.add_argument("input_pod5s", nargs="+", type=Path, help="Input .pod5 files")
    return ap

def main(argv: Optional[List[str]] = None) -> int:
    args = build_parser().parse_args(argv)
    if args.output_pod5.exists():
        raise FileExistsError(f"Output file already exists: {args.output_pod5}")
    if args.output_pod5.suffix.lower() != ".pod5":
        raise ValueError(f"Output does not look like a .pod5 file: {args.output_pod5}")

    for input_pod5 in args.input_pod5s:
        if not input_pod5.exists():
            raise FileNotFoundError(f"Input file does not exist: {input_pod5}")
        if input_pod5.suffix.lower() != ".pod5":
            raise ValueError(f"Input does not look like a .pod5 file: {input_pod5}")

    counter = 0;
    readers = [p5.Reader(str(pod5_path)) for pod5_path in args.input_pod5s]
    with p5.Writer(str(args.output_pod5)) as writer:
        # Merge read record from all input POD5 files by taken 1 from each in round-robin fashion
        reader_iterators = [iter(reader.reads()) for reader in readers]
        while any(reader_iterators):
            for reader_iterator in reader_iterators:
                try:
                    read_record = next(reader_iterator)
                    writer.add_read(read_record.to_read())
                    counter += 1
                except StopIteration:
                    reader_iterators.remove(reader_iterator)
    for reader in readers:
        reader.close()
    print(f"Merged {counter} reads into {args.output_pod5}")
    return 0

if __name__ == "__main__":
    raise SystemExit(main())
