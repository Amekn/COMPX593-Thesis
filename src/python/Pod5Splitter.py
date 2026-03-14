#!/usr/bin/env python3
"""Split a POD5 file into multiple outputs with balanced read counts."""

from __future__ import annotations

import argparse
from contextlib import ExitStack
from pathlib import Path
from typing import Sequence

import pod5 as p5


def build_argument_parser() -> argparse.ArgumentParser:
    """Create the command-line interface for the POD5 splitter."""
    parser = argparse.ArgumentParser(
        description="Split one POD5 file into multiple POD5 outputs with near-equal read counts."
    )
    parser.add_argument("source_pod5", type=Path, help="Input POD5 file.")
    parser.add_argument("number_of_output_files", type=int, help="Number of output POD5 files.")
    parser.add_argument("output_directory", type=Path, help="Destination directory.")
    parser.add_argument("output_prefix", help="Prefix to use for each output filename.")
    return parser


def build_output_paths(output_directory: Path, output_prefix: str, output_count: int) -> list[Path]:
    """Construct deterministic output filenames for all split files."""
    return [
        output_directory / f"{output_prefix}.part{output_index + 1:04d}.pod5"
        for output_index in range(output_count)
    ]


def compute_partition_sizes(total_reads: int, output_count: int) -> list[int]:
    """Distribute reads as evenly as possible across the requested number of outputs."""
    base_size, remainder = divmod(total_reads, output_count)
    return [
        base_size + (1 if output_index < remainder else 0)
        for output_index in range(output_count)
    ]


def validate_arguments(source_path: Path, output_paths: Sequence[Path], output_count: int) -> None:
    """Validate the requested split before any output files are created."""
    if not source_path.exists():
        raise FileNotFoundError(f"Input file not found: {source_path}")
    if source_path.suffix.lower() != ".pod5":
        raise ValueError(f"Input path must use the .pod5 extension: {source_path}")
    if output_count <= 0:
        raise ValueError("number_of_output_files must be at least 1")
    if any(output_path.exists() for output_path in output_paths):
        existing_paths = ", ".join(str(output_path) for output_path in output_paths if output_path.exists())
        raise FileExistsError(f"Refusing to overwrite existing output files: {existing_paths}")


def split_pod5(source_path: Path, output_paths: Sequence[Path]) -> list[int]:
    """Write balanced POD5 splits and return the number of reads written to each output file."""
    with ExitStack() as stack:
        reader = stack.enter_context(p5.Reader(str(source_path)))
        total_reads = int(reader.num_reads)
        if total_reads == 0:
            raise RuntimeError("Input POD5 contains no reads.")
        if len(output_paths) > total_reads:
            raise ValueError(
                f"Requested {len(output_paths)} outputs for {total_reads} reads. "
                "This would create empty output files."
            )

        target_sizes = compute_partition_sizes(total_reads, len(output_paths))
        writers = [stack.enter_context(p5.Writer(str(output_path))) for output_path in output_paths]
        written_counts = [0 for _ in output_paths]

        current_writer_index = 0
        reads_written_to_current_output = 0

        for read_record in reader.reads():
            writers[current_writer_index].add_read(read_record.to_read())
            written_counts[current_writer_index] += 1
            reads_written_to_current_output += 1

            if reads_written_to_current_output == target_sizes[current_writer_index]:
                current_writer_index += 1
                reads_written_to_current_output = 0
                if current_writer_index == len(writers):
                    break

    return written_counts


def main(argv: Sequence[str] | None = None) -> int:
    """Run the command-line entry point."""
    arguments = build_argument_parser().parse_args(argv)

    source_path = arguments.source_pod5.resolve()
    output_directory = arguments.output_directory.resolve()
    output_paths = build_output_paths(
        output_directory=output_directory,
        output_prefix=arguments.output_prefix,
        output_count=arguments.number_of_output_files,
    )

    validate_arguments(source_path, output_paths, arguments.number_of_output_files)
    output_directory.mkdir(parents=True, exist_ok=True)

    written_counts = split_pod5(source_path, output_paths)
    for output_path, written_count in zip(output_paths, written_counts):
        print(f"{output_path}\t{written_count}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
