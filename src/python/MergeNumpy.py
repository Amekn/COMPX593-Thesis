#!/usr/bin/env python3
"""Merge Bonito training arrays from multiple directories without large RAM spikes."""

from __future__ import annotations

import argparse
import math
import os
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable, Sequence

import numpy as np
from numpy.lib.format import open_memmap

try:
    from tqdm import tqdm as progress_bar
except Exception:  # pragma: no cover
    class _NullProgressBar:
        """Fallback progress bar used when tqdm is unavailable."""

        def __init__(self, *_args, **_kwargs) -> None:
            pass

        def update(self, _increment: int = 1) -> None:
            pass

        def close(self) -> None:
            pass

    def progress_bar(*args, **kwargs):  # type: ignore[misc]
        return _NullProgressBar(*args, **kwargs)


REQUIRED_FILENAMES = ("chunks.npy", "references.npy", "reference_lengths.npy")


class MergeError(RuntimeError):
    """Raised when the input datasets cannot be merged safely."""


@dataclass(frozen=True)
class ArrayTripletPaths:
    """Locations of the three required training arrays inside one dataset directory."""

    chunks: Path
    references: Path
    reference_lengths: Path


@dataclass(frozen=True)
class DatasetMetadata:
    """Lightweight header information gathered from one input dataset."""

    directory: Path
    row_count: int
    chunk_tail_shape: tuple[int, ...]
    reference_width: int
    chunks_dtype: np.dtype
    references_dtype: np.dtype
    reference_lengths_dtype: np.dtype


@dataclass(frozen=True)
class OutputPaths:
    """Destination filenames for the merged training arrays."""

    chunks: Path
    references: Path
    reference_lengths: Path


@dataclass(frozen=True)
class MergePlan:
    """Derived merge layout shared across the streaming copy and shuffle stages."""

    datasets: list[DatasetMetadata]
    total_rows: int
    chunk_shape: tuple[int, ...]
    reference_shape: tuple[int, int]
    reference_length_shape: tuple[int]
    chunks_dtype: np.dtype
    references_dtype: np.dtype
    reference_lengths_dtype: np.dtype


def build_argument_parser() -> argparse.ArgumentParser:
    """Create the command-line interface for the array merger."""
    parser = argparse.ArgumentParser(
        description=(
            "Merge Bonito training arrays (chunks, references, reference_lengths) "
            "from multiple directories without loading every array into RAM."
        )
    )
    parser.add_argument(
        "input_dirs",
        nargs="+",
        help="Input directories. Each must contain chunks.npy, references.npy, and reference_lengths.npy.",
    )
    parser.add_argument("-o", "--output", required=True, help="Output directory.")
    parser.add_argument(
        "--overwrite",
        action="store_true",
        help="Allow replacement of existing output array files.",
    )
    parser.add_argument(
        "--target-chunk-mb",
        type=float,
        default=64.0,
        help="Approximate RAM budget, in MiB, for each streaming block.",
    )
    parser.add_argument(
        "--shuffle",
        action="store_true",
        help="Shuffle rows after merging while preserving row correspondence.",
    )
    parser.add_argument(
        "--seed",
        type=int,
        default=None,
        help="Random seed used when --shuffle is enabled.",
    )
    return parser


def resolve_required_paths(dataset_directory: Path) -> ArrayTripletPaths:
    """Locate the three required NumPy files in one input directory."""
    paths = ArrayTripletPaths(
        chunks=dataset_directory / REQUIRED_FILENAMES[0],
        references=dataset_directory / REQUIRED_FILENAMES[1],
        reference_lengths=dataset_directory / REQUIRED_FILENAMES[2],
    )

    missing_files = [
        path.name
        for path in (paths.chunks, paths.references, paths.reference_lengths)
        if not path.exists()
    ]
    if missing_files:
        missing_file_list = ", ".join(missing_files)
        raise MergeError(f"Missing required files in {dataset_directory}: {missing_file_list}")

    return paths


def inspect_dataset(dataset_directory: Path) -> DatasetMetadata:
    """Read array headers with memmap mode and return merge-relevant metadata."""
    triplet_paths = resolve_required_paths(dataset_directory)
    chunk_array = np.load(triplet_paths.chunks, allow_pickle=False, mmap_mode="r")
    reference_array = np.load(triplet_paths.references, allow_pickle=False, mmap_mode="r")
    reference_length_array = np.load(
        triplet_paths.reference_lengths, allow_pickle=False, mmap_mode="r"
    )

    if chunk_array.ndim < 1:
        raise MergeError(f"chunks.npy must be at least 1-D, got shape {chunk_array.shape}")
    if reference_array.ndim != 2:
        raise MergeError(f"references.npy must be 2-D, got shape {reference_array.shape}")
    if reference_length_array.ndim != 1:
        raise MergeError(
            "reference_lengths.npy must be 1-D, "
            f"got shape {reference_length_array.shape}"
        )

    row_count = int(chunk_array.shape[0])
    if row_count != int(reference_array.shape[0]) or row_count != int(reference_length_array.shape[0]):
        raise MergeError(
            "Row-count mismatch in "
            f"{dataset_directory} (chunks={chunk_array.shape[0]}, "
            f"references={reference_array.shape[0]}, "
            f"reference_lengths={reference_length_array.shape[0]})"
        )

    return DatasetMetadata(
        directory=dataset_directory,
        row_count=row_count,
        chunk_tail_shape=tuple(int(dimension) for dimension in chunk_array.shape[1:]),
        reference_width=int(reference_array.shape[1]),
        chunks_dtype=np.dtype(chunk_array.dtype),
        references_dtype=np.dtype(reference_array.dtype),
        reference_lengths_dtype=np.dtype(reference_length_array.dtype),
    )


def compute_common_dtype(dtypes: Sequence[np.dtype]) -> np.dtype:
    """Return a common dtype that can safely represent all supplied arrays."""
    common = dtypes[0]
    for dtype in dtypes[1:]:
        common = np.result_type(common, dtype)
    return np.dtype(common)


def build_merge_plan(input_directories: Sequence[Path]) -> MergePlan:
    """Inspect all inputs and derive the merged output shape and dtypes."""
    dataset_metadata: list[DatasetMetadata] = []
    progress = progress_bar(total=len(input_directories), desc="Inspecting headers", unit="dir")

    try:
        for dataset_directory in input_directories:
            if not dataset_directory.is_dir():
                raise MergeError(f"Input path is not a directory: {dataset_directory}")
            dataset_metadata.append(inspect_dataset(dataset_directory))
            progress.update(1)
    finally:
        progress.close()

    if not dataset_metadata:
        raise MergeError("At least one input directory is required.")

    chunk_tail_shapes = {metadata.chunk_tail_shape for metadata in dataset_metadata}
    if len(chunk_tail_shapes) != 1:
        shape_list = ", ".join(str(shape) for shape in sorted(chunk_tail_shapes))
        raise MergeError(f"Inconsistent chunk tail shapes across inputs: {shape_list}")

    total_rows = sum(metadata.row_count for metadata in dataset_metadata)
    max_reference_width = max(metadata.reference_width for metadata in dataset_metadata)
    chunk_tail_shape = next(iter(chunk_tail_shapes))

    chunks_dtype = compute_common_dtype([metadata.chunks_dtype for metadata in dataset_metadata])
    references_dtype = compute_common_dtype(
        [metadata.references_dtype for metadata in dataset_metadata]
    )
    reference_lengths_dtype = compute_common_dtype(
        [metadata.reference_lengths_dtype for metadata in dataset_metadata]
    )

    return MergePlan(
        datasets=dataset_metadata,
        total_rows=total_rows,
        chunk_shape=(total_rows, *chunk_tail_shape),
        reference_shape=(total_rows, max_reference_width),
        reference_length_shape=(total_rows,),
        chunks_dtype=chunks_dtype,
        references_dtype=references_dtype,
        reference_lengths_dtype=reference_lengths_dtype,
    )


def resolve_output_paths(output_directory: Path) -> OutputPaths:
    """Construct the destination array paths in the output directory."""
    return OutputPaths(
        chunks=output_directory / REQUIRED_FILENAMES[0],
        references=output_directory / REQUIRED_FILENAMES[1],
        reference_lengths=output_directory / REQUIRED_FILENAMES[2],
    )


def ensure_output_paths_are_writable(output_paths: OutputPaths, overwrite: bool) -> None:
    """Prevent accidental replacement of an existing merged dataset unless requested."""
    for output_path in (output_paths.chunks, output_paths.references, output_paths.reference_lengths):
        if output_path.exists() and not overwrite:
            raise MergeError(
                f"Output file already exists: {output_path}. Use --overwrite to replace it."
            )


def create_output_memmaps(plan: MergePlan, output_paths: OutputPaths) -> tuple[np.memmap, np.memmap, np.memmap]:
    """Allocate the merged arrays as writable memmaps on disk."""
    print("[INFO] Creating output memmaps ...")
    merged_chunks = open_memmap(
        output_paths.chunks,
        mode="w+",
        dtype=plan.chunks_dtype,
        shape=plan.chunk_shape,
    )
    merged_references = open_memmap(
        output_paths.references,
        mode="w+",
        dtype=plan.references_dtype,
        shape=plan.reference_shape,
    )
    merged_reference_lengths = open_memmap(
        output_paths.reference_lengths,
        mode="w+",
        dtype=plan.reference_lengths_dtype,
        shape=plan.reference_length_shape,
    )

    # Reference rows are zero-padded out to the global maximum width.
    merged_references[:] = 0
    return merged_chunks, merged_references, merged_reference_lengths


def estimate_rows_per_block(
    metadata: DatasetMetadata,
    plan: MergePlan,
    target_bytes: int,
) -> int:
    """Choose a row block size that fits approximately within the target memory budget."""
    chunk_elements_per_row = int(math.prod(metadata.chunk_tail_shape)) if metadata.chunk_tail_shape else 1
    chunk_bytes_per_row = chunk_elements_per_row * int(np.dtype(plan.chunks_dtype).itemsize)
    reference_bytes_per_row = metadata.reference_width * int(np.dtype(plan.references_dtype).itemsize)
    reference_length_bytes_per_row = int(np.dtype(plan.reference_lengths_dtype).itemsize)

    bytes_per_row = chunk_bytes_per_row + reference_bytes_per_row + reference_length_bytes_per_row
    return max(1, min(metadata.row_count, target_bytes // max(1, bytes_per_row)))


def stream_merge(
    plan: MergePlan,
    merged_chunks: np.memmap,
    merged_references: np.memmap,
    merged_reference_lengths: np.memmap,
    target_chunk_mb: float,
) -> None:
    """Copy all input rows into the output memmaps in bounded blocks."""
    print("[INFO] Streaming arrays into output (mmap -> mmap) ...")
    target_bytes = max(1, int(target_chunk_mb * 1024 * 1024))
    progress = progress_bar(total=plan.total_rows, desc="Streaming rows", unit="row")

    destination_row_start = 0
    try:
        for metadata in plan.datasets:
            triplet_paths = resolve_required_paths(metadata.directory)
            source_chunks = np.load(triplet_paths.chunks, allow_pickle=False, mmap_mode="r")
            source_references = np.load(triplet_paths.references, allow_pickle=False, mmap_mode="r")
            source_reference_lengths = np.load(
                triplet_paths.reference_lengths, allow_pickle=False, mmap_mode="r"
            )

            rows_per_block = estimate_rows_per_block(metadata, plan, target_bytes)
            source_row_start = 0
            while source_row_start < metadata.row_count:
                source_row_end = min(metadata.row_count, source_row_start + rows_per_block)
                destination_row_end = destination_row_start + (source_row_end - source_row_start)

                merged_chunks[destination_row_start:destination_row_end, ...] = source_chunks[
                    source_row_start:source_row_end, ...
                ]
                merged_references[
                    destination_row_start:destination_row_end, : metadata.reference_width
                ] = source_references[source_row_start:source_row_end, : metadata.reference_width]
                merged_reference_lengths[destination_row_start:destination_row_end] = (
                    source_reference_lengths[source_row_start:source_row_end]
                )

                progress.update(source_row_end - source_row_start)
                destination_row_start = destination_row_end
                source_row_start = source_row_end
    finally:
        progress.close()

    merged_chunks.flush()
    merged_references.flush()
    merged_reference_lengths.flush()


def shuffle_output_rows(
    plan: MergePlan,
    output_paths: OutputPaths,
    target_chunk_mb: float,
    seed: int | None,
) -> None:
    """Shuffle the merged arrays on disk while keeping row alignment intact."""
    print("[INFO] Shuffling merged rows ...")
    merged_chunks = np.load(output_paths.chunks, allow_pickle=False, mmap_mode="r")
    merged_references = np.load(output_paths.references, allow_pickle=False, mmap_mode="r")
    merged_reference_lengths = np.load(output_paths.reference_lengths, allow_pickle=False, mmap_mode="r")

    rng = np.random.default_rng(seed)
    permutation = rng.permutation(plan.total_rows)

    temporary_paths = OutputPaths(
        chunks=output_paths.chunks.with_name("chunks.tmp.npy"),
        references=output_paths.references.with_name("references.tmp.npy"),
        reference_lengths=output_paths.reference_lengths.with_name("reference_lengths.tmp.npy"),
    )

    shuffled_chunks = open_memmap(
        temporary_paths.chunks,
        mode="w+",
        dtype=plan.chunks_dtype,
        shape=plan.chunk_shape,
    )
    shuffled_references = open_memmap(
        temporary_paths.references,
        mode="w+",
        dtype=plan.references_dtype,
        shape=plan.reference_shape,
    )
    shuffled_reference_lengths = open_memmap(
        temporary_paths.reference_lengths,
        mode="w+",
        dtype=plan.reference_lengths_dtype,
        shape=plan.reference_length_shape,
    )

    chunk_elements_per_row = int(math.prod(plan.chunk_shape[1:])) if len(plan.chunk_shape) > 1 else 1
    chunk_bytes_per_row = chunk_elements_per_row * int(np.dtype(plan.chunks_dtype).itemsize)
    reference_bytes_per_row = plan.reference_shape[1] * int(np.dtype(plan.references_dtype).itemsize)
    reference_length_bytes_per_row = int(np.dtype(plan.reference_lengths_dtype).itemsize)
    bytes_per_row = chunk_bytes_per_row + reference_bytes_per_row + reference_length_bytes_per_row
    rows_per_block = max(1, int(max(1, target_chunk_mb * 1024 * 1024) // max(1, bytes_per_row)))

    progress = progress_bar(total=plan.total_rows, desc="Shuffling", unit="row")
    try:
        for shuffled_row_start in range(0, plan.total_rows, rows_per_block):
            shuffled_row_end = min(plan.total_rows, shuffled_row_start + rows_per_block)
            row_indices = permutation[shuffled_row_start:shuffled_row_end]

            shuffled_chunks[shuffled_row_start:shuffled_row_end, ...] = merged_chunks[row_indices, ...]
            shuffled_references[shuffled_row_start:shuffled_row_end, ...] = (
                merged_references[row_indices, ...]
            )
            shuffled_reference_lengths[shuffled_row_start:shuffled_row_end] = (
                merged_reference_lengths[row_indices]
            )
            progress.update(shuffled_row_end - shuffled_row_start)
    finally:
        progress.close()

    shuffled_chunks.flush()
    shuffled_references.flush()
    shuffled_reference_lengths.flush()

    os.replace(temporary_paths.chunks, output_paths.chunks)
    os.replace(temporary_paths.references, output_paths.references)
    os.replace(temporary_paths.reference_lengths, output_paths.reference_lengths)


def print_summary(plan: MergePlan, output_paths: OutputPaths, shuffled: bool) -> None:
    """Print a concise summary of the merged dataset written to disk."""
    print(f"[DONE] Wrote merged arrays to {output_paths.chunks.parent}:")
    if shuffled:
        print("       rows shuffled: yes")
    print(
        f"       {output_paths.chunks.name}: shape {plan.chunk_shape}, "
        f"dtype {plan.chunks_dtype}"
    )
    print(
        f"       {output_paths.references.name}: shape {plan.reference_shape}, "
        f"dtype {plan.references_dtype}"
    )
    print(
        f"       {output_paths.reference_lengths.name}: shape {plan.reference_length_shape}, "
        f"dtype {plan.reference_lengths_dtype}"
    )
    print(f"[DONE] Total records: {plan.total_rows}")


def main(argv: Sequence[str] | None = None) -> int:
    """Run the command-line entry point."""
    arguments = build_argument_parser().parse_args(argv)
    input_directories = [Path(directory).resolve() for directory in arguments.input_dirs]
    output_directory = Path(arguments.output).resolve()
    output_directory.mkdir(parents=True, exist_ok=True)

    print(
        f"[INFO] Inspecting {len(input_directories)} input director"
        f"{'y' if len(input_directories) == 1 else 'ies'} ..."
    )
    plan = build_merge_plan(input_directories)
    output_paths = resolve_output_paths(output_directory)
    ensure_output_paths_are_writable(output_paths, overwrite=arguments.overwrite)

    merged_chunks, merged_references, merged_reference_lengths = create_output_memmaps(
        plan=plan,
        output_paths=output_paths,
    )
    stream_merge(
        plan=plan,
        merged_chunks=merged_chunks,
        merged_references=merged_references,
        merged_reference_lengths=merged_reference_lengths,
        target_chunk_mb=arguments.target_chunk_mb,
    )

    if arguments.shuffle:
        shuffle_output_rows(
            plan=plan,
            output_paths=output_paths,
            target_chunk_mb=arguments.target_chunk_mb,
            seed=arguments.seed,
        )

    print_summary(plan, output_paths, shuffled=arguments.shuffle)
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except MergeError as exception:
        print(f"[ERROR] {exception}", file=sys.stderr)
        raise SystemExit(2)
