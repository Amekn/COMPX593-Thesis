#!/usr/bin/env python3
"""
pod5Spliter.py

Split a single POD5 file into multiple POD5 files with (nearly) equal numbers of reads.

Requires:
  pip install pod5

Example:
  python pod5Spliter.py <source_pod5> <number_of_output_files> <output_directory> <output_prefix>
"""

from __future__ import annotations

import argparse
from pathlib import Path
from typing import List, Optional

import pod5 as p5

def build_parser() -> argparse.ArgumentParser:
    ap = argparse.ArgumentParser(
        description="Split a POD5 file into multiple outputs with (nearly) equal reads per file."
    )
    ap.add_argument("source_pod5", type=Path, help="Input .pod5 file")
    ap.add_argument("number_of_output_files", type=int, help="Number of output .pod5 files to create")
    ap.add_argument("output_directory", type=Path, help="Output directory")
    ap.add_argument("output_prefix", type=str, help="Output filename prefix")
    return ap

def main(argv: Optional[List[str]] = None) -> int:
    args = build_parser().parse_args(argv)

    if not args.source_pod5.exists():
        raise FileNotFoundError(f"Input file not found: {args.source_pod5}")
    if args.source_pod5.suffix.lower() != ".pod5":
        raise ValueError(f"Input does not look like a .pod5 file: {args.source_pod5}")

    if args.number_of_output_files <= 0:
        raise ValueError("number_of_output_files must be >= 1")

    outdir: Path = args.output_directory
    outdir.mkdir(parents=True, exist_ok=True)

    prefix = args.output_prefix

    writers = []
    for i in range(args.number_of_output_files):
        out_path = outdir / f"{prefix}.part{i+1:04d}.pod5"
        writers.append(p5.Writer(str(out_path)))
        
    with p5.Reader(str(args.source_pod5)) as reader:
        total_reads = reader.num_reads

        if total_reads == 0:
            raise RuntimeError("Input POD5 contains zero reads. Congrats on your very empty file.")

        if args.number_of_output_files > total_reads:
            raise ValueError(
                f"--number_of_output_files ({args.number_of_output_files}) is greater than total reads ({total_reads}). "
                "That would create empty outputs, which is not 'equal' in any meaningful sense."
            )
        
        batch_size = total_reads / args.number_of_output_files

        counter = 0
        current_writer_index = 0
        for read_record in reader.reads():
            if counter <= batch_size:
                writers[current_writer_index].add_read(read_record.to_read())
                counter += 1
            else: 
                print(f"Finished writing {counter} reads to {writers[current_writer_index].path}")
                writers[current_writer_index].close()
                current_writer_index += 1
                counter = 0
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
