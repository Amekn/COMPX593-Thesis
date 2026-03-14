#!/usr/bin/env python3
"""Merge multiple POD5 files into a single output file."""

from __future__ import annotations

import argparse
from contextlib import ExitStack
from pathlib import Path
from typing import Sequence

import pod5 as p5


def build_argument_parser() -> argparse.ArgumentParser:
    """Create the command-line interface for the POD5 merger."""
    parser = argparse.ArgumentParser(
        description=(
            "Merge multiple POD5 inputs into a single output file. "
            "Reads are written in round-robin order so each source contributes evenly."
        )
    )
    parser.add_argument("output_pod5", type=Path, help="Destination POD5 file.")
    parser.add_argument("input_pod5s", nargs="+", type=Path, help="Input POD5 files to merge.")
    return parser


def validate_paths(output_path: Path, input_paths: Sequence[Path]) -> None:
    """Validate that the requested merge paths are sensible before writing anything."""
    if output_path.exists():
        raise FileExistsError(f"Output file already exists: {output_path}")
    if output_path.suffix.lower() != ".pod5":
        raise ValueError(f"Output path must use the .pod5 extension: {output_path}")

    resolved_output = output_path.resolve()
    resolved_inputs = []
    for input_path in input_paths:
        if not input_path.exists():
            raise FileNotFoundError(f"Input file does not exist: {input_path}")
        if input_path.suffix.lower() != ".pod5":
            raise ValueError(f"Input path must use the .pod5 extension: {input_path}")
        resolved_inputs.append(input_path.resolve())

    if resolved_output in resolved_inputs:
        raise ValueError("Output POD5 path must not match any input POD5 path.")


def merge_round_robin(output_path: Path, input_paths: Sequence[Path]) -> int:
    """Write reads from each input file in round-robin order and return the count written."""
    with ExitStack() as stack:
        readers = [stack.enter_context(p5.Reader(str(input_path))) for input_path in input_paths]
        writer = stack.enter_context(p5.Writer(str(output_path)))

        active_iterators = [iter(reader.reads()) for reader in readers]
        written_read_count = 0

        while active_iterators:
            next_round_iterators = []
            for read_iterator in active_iterators:
                try:
                    read_record = next(read_iterator)
                except StopIteration:
                    continue

                writer.add_read(read_record.to_read())
                written_read_count += 1
                next_round_iterators.append(read_iterator)

            active_iterators = next_round_iterators

    return written_read_count


def main(argv: Sequence[str] | None = None) -> int:
    """Run the command-line entry point."""
    arguments = build_argument_parser().parse_args(argv)
    output_path = arguments.output_pod5.resolve()
    input_paths = [input_path.resolve() for input_path in arguments.input_pod5s]

    validate_paths(output_path, input_paths)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    written_read_count = merge_round_robin(output_path, input_paths)
    print(f"Merged {written_read_count} reads into {output_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
