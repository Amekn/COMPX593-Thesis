#!/usr/bin/env python3
"""Merge multiple POD5 files into a single POD5 output using round-robin ordering.

Reassembles POD5 shards produced by ``Pod5Splitter.py`` or combines separate
A/B datasets in the COMPX593 thesis pipeline. Rather than concatenating each
input end-to-end, reads are interleaved one-at-a-time across all inputs so
the resulting file draws evenly from every source even when downstream
consumers only process a prefix. The merger never rewrites metadata; read
records are copied verbatim via ``read_record.to_read()``.
"""

from __future__ import annotations

import argparse
from contextlib import ExitStack
from pathlib import Path
from typing import Sequence

import pod5 as p5


def build_argument_parser() -> argparse.ArgumentParser:
    """Create the CLI parser for the POD5 merger.

    The first positional argument is the destination path; all subsequent
    positional arguments are POD5 inputs. At least one input is required.
    """
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
    """Check merge paths are sane before any output file is opened.

    Rejects the request if the output already exists, if any path has the
    wrong extension, if an input is missing, or if the resolved output path
    also appears in the input list (which would truncate and corrupt the
    source mid-merge).

    Raises:
        FileExistsError: If the output path already exists.
        FileNotFoundError: If any input path is missing.
        ValueError: If extensions are wrong or an input aliases the output.
    """
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
    """Interleave reads from every input into ``output_path`` round-robin.

    On each pass the function pulls exactly one read from every input whose
    iterator has not yet been exhausted. Exhausted iterators drop out of the
    active set for the next pass, so the final tail of the merged file is
    simply the leftover reads of whichever inputs were longest.

    Args:
        output_path: Destination POD5 path. Must not already exist.
        input_paths: One or more POD5 inputs to draw from.

    Returns:
        The total number of reads written.
    """
    # ExitStack guarantees every reader and the writer are closed even if a
    # mid-stream failure raises, so partial output files still flush cleanly.
    with ExitStack() as stack:
        readers = [stack.enter_context(p5.Reader(str(input_path))) for input_path in input_paths]
        writer = stack.enter_context(p5.Writer(str(output_path)))

        active_iterators = [iter(reader.reads()) for reader in readers]
        written_read_count = 0

        while active_iterators:
            # Rebuild the active set each pass; inputs that raised
            # StopIteration are silently dropped for subsequent rounds.
            next_round_iterators = []
            for read_iterator in active_iterators:
                try:
                    read_record = next(read_iterator)
                except StopIteration:
                    continue

                # ``to_read()`` materialises an owned, mutable copy; required
                # because the pod5 Writer rejects streaming ReadRecord views.
                writer.add_read(read_record.to_read())
                written_read_count += 1
                next_round_iterators.append(read_iterator)

            active_iterators = next_round_iterators

    return written_read_count


def main(argv: Sequence[str] | None = None) -> int:
    """CLI entry point: validate paths, run the merge, and report the total.

    Args:
        argv: Optional argument vector for testing; defaults to ``sys.argv[1:]``.

    Returns:
        Process exit status (always ``0`` on success).
    """
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
