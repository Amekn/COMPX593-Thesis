#!/usr/bin/env python3
"""Partition a POD5 file into N shards with near-equal read counts.

Used in the COMPX593 thesis preprocessing stage to shard a large raw POD5
so each shard can be basecalled in parallel on a separate GPU/worker. Reads
are emitted in source order: the first ``k`` reads fill shard 0, the next
``k`` fill shard 1, and so on, where shard sizes differ by at most one read
when ``total_reads`` is not evenly divisible. Read records are copied
verbatim; no metadata is rewritten.
"""

from __future__ import annotations

import argparse
from contextlib import ExitStack
from pathlib import Path
from typing import Sequence

import pod5 as p5


def build_argument_parser() -> argparse.ArgumentParser:
    """Create the CLI parser for the POD5 splitter.

    Accepts the input POD5, the requested shard count, the output directory,
    and a filename prefix used to build deterministic shard names.
    """
    parser = argparse.ArgumentParser(
        description="Split one POD5 file into multiple POD5 outputs with near-equal read counts."
    )
    parser.add_argument("source_pod5", type=Path, help="Input POD5 file.")
    parser.add_argument("number_of_output_files", type=int, help="Number of output POD5 files.")
    parser.add_argument("output_directory", type=Path, help="Destination directory.")
    parser.add_argument("output_prefix", help="Prefix to use for each output filename.")
    return parser


def build_output_paths(output_directory: Path, output_prefix: str, output_count: int) -> list[Path]:
    """Construct deterministic shard filenames of the form ``<prefix>.partNNNN.pod5``.

    The four-digit ``%04d`` index provides lexicographic ordering that
    matches numeric ordering for up to 9999 shards, which is far beyond any
    realistic parallel-basecall batch.
    """
    return [
        output_directory / f"{output_prefix}.part{output_index + 1:04d}.pod5"
        for output_index in range(output_count)
    ]


def compute_partition_sizes(total_reads: int, output_count: int) -> list[int]:
    """Return shard sizes summing to ``total_reads`` and differing by at most one.

    The first ``remainder`` shards receive an extra read, matching the
    classic "ceil for the first r, floor for the rest" balanced-split rule.
    """
    base_size, remainder = divmod(total_reads, output_count)
    return [
        base_size + (1 if output_index < remainder else 0)
        for output_index in range(output_count)
    ]


def validate_arguments(source_path: Path, output_paths: Sequence[Path], output_count: int) -> None:
    """Check the requested split before any shard file is opened for writing.

    Raises:
        FileNotFoundError: If the source POD5 is missing.
        ValueError: If the source extension is wrong or ``output_count <= 0``.
        FileExistsError: If any proposed shard path already exists.
    """
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
    """Stream the source POD5 into the prepared shard writers.

    Opens all writers up front through a single ``ExitStack`` so every shard
    is flushed and closed even when a mid-stream error occurs. Reads are
    distributed contiguously: shard ``i`` receives the ``i``-th contiguous
    block of ``target_sizes[i]`` reads in source-file order.

    Args:
        source_path: POD5 file to partition.
        output_paths: Pre-computed shard destination paths.

    Returns:
        List of per-shard read counts in the same order as ``output_paths``.

    Raises:
        RuntimeError: If the source POD5 has no reads.
        ValueError: If more shards were requested than there are reads.
    """
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

        # Rolling cursor: we advance ``current_writer_index`` only when the
        # current shard has received exactly its target size, which keeps
        # reads contiguous within each shard (simplifies cross-referencing
        # against the original file when debugging a specific basecall).
        current_writer_index = 0
        reads_written_to_current_output = 0

        for read_record in reader.reads():
            writers[current_writer_index].add_read(read_record.to_read())
            written_counts[current_writer_index] += 1
            reads_written_to_current_output += 1

            if reads_written_to_current_output == target_sizes[current_writer_index]:
                current_writer_index += 1
                reads_written_to_current_output = 0
                # Break early: any reads remaining in the source beyond the
                # planned total would overflow past the last shard.
                if current_writer_index == len(writers):
                    break

    return written_counts


def main(argv: Sequence[str] | None = None) -> int:
    """CLI entry point: resolve paths, validate, split, and print counts.

    On success prints one ``<path>\\t<count>`` line per shard to stdout so
    the tab-separated output composes with shell pipelines.
    """
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
