#!/usr/bin/env python3
"""Memory-map merge of Bonito training-array triplets without RAM spikes.

Bonito exports each training split as three parallel NumPy files --
``chunks.npy``, ``references.npy``, ``reference_lengths.npy`` -- whose
leading axis is a shared row (chunk) count. This utility concatenates
multiple such triplets on disk using ``numpy.lib.format.open_memmap`` so
building a combined A+B training set for the COMPX593 thesis fine-tune
never has to hold all inputs in RAM at once.

Key invariants the implementation enforces:
    * Every input dataset's three arrays share the same leading axis.
    * The ``chunks`` tail shape (all axes after the leading one) is
      identical across every input; mismatches abort the merge.
    * ``references`` widths may differ per input; the merged array is
      padded to the global maximum width with zeros.
    * Output dtypes are promoted via ``np.result_type`` so no input value
      is silently narrowed.

Optional features:
    * ``--shuffle`` performs an on-disk permutation after the streaming
      concat, keeping the three arrays' row correspondence intact.
    * ``--target-chunk-mb`` bounds the per-block RAM budget of every
      mmap-to-mmap copy.
"""

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
    # tqdm is optional; fall back to a silent shim so the tool runs on
    # minimal training-box environments without pulling in extras.
    class _NullProgressBar:
        """No-op stand-in for tqdm used when the library is not installed."""

        def __init__(self, *_args, **_kwargs) -> None:
            pass

        def update(self, _increment: int = 1) -> None:
            pass

        def close(self) -> None:
            pass

    def progress_bar(*args, **kwargs):  # type: ignore[misc]
        return _NullProgressBar(*args, **kwargs)


# Bonito's own training loader hard-codes these three filenames, so every
# input dataset must expose them verbatim and every output directory
# produces the same trio to remain drop-in compatible.
REQUIRED_FILENAMES = ("chunks.npy", "references.npy", "reference_lengths.npy")


class MergeError(RuntimeError):
    """Raised when the inputs fail validation or the merge cannot proceed.

    Distinguished from unrelated ``RuntimeError``s so the top-level
    ``__main__`` block can print a concise ``[ERROR]`` line and exit with
    status 2 rather than dumping a traceback.
    """


@dataclass(frozen=True)
class ArrayTripletPaths:
    """Resolved paths of the three required arrays inside one input directory."""

    chunks: Path
    references: Path
    reference_lengths: Path


@dataclass(frozen=True)
class DatasetMetadata:
    """Per-input header snapshot used to plan the merge without reading payloads.

    Stores row count, the invariant ``chunks`` tail shape, the per-dataset
    reference width (may differ across inputs), and each array's dtype.
    """

    directory: Path
    row_count: int
    chunk_tail_shape: tuple[int, ...]
    reference_width: int
    chunks_dtype: np.dtype
    references_dtype: np.dtype
    reference_lengths_dtype: np.dtype


@dataclass(frozen=True)
class OutputPaths:
    """Resolved destination paths for the three merged arrays."""

    chunks: Path
    references: Path
    reference_lengths: Path


@dataclass(frozen=True)
class MergePlan:
    """Computed merge layout: totals, output shapes, and promoted dtypes.

    Produced once up front by ``build_merge_plan`` and then consumed by the
    streaming-copy and optional shuffle stages, guaranteeing they agree on
    total row count, output shapes, and dtype promotions.
    """

    datasets: list[DatasetMetadata]
    total_rows: int
    chunk_shape: tuple[int, ...]
    reference_shape: tuple[int, int]
    reference_length_shape: tuple[int]
    chunks_dtype: np.dtype
    references_dtype: np.dtype
    reference_lengths_dtype: np.dtype


def build_argument_parser() -> argparse.ArgumentParser:
    """Create the CLI parser for the training-array merger.

    At least one input directory is required. ``--target-chunk-mb`` caps
    the RAM footprint of each row-block copy; 64 MiB is a conservative
    default that streams comfortably on training nodes with many other
    processes active. ``--shuffle`` runs an extra on-disk permutation pass
    (row-aligned across the three arrays) and ``--seed`` makes that pass
    deterministic for reproducibility.
    """
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
    """Return the triplet paths for a dataset, raising if any file is missing.

    Raises:
        MergeError: If one or more of ``chunks.npy``, ``references.npy``,
            or ``reference_lengths.npy`` is absent.
    """
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
    """Collect header-only metadata for one dataset without reading payloads.

    Uses ``mmap_mode="r"`` so only the ``.npy`` header plus a mapped view
    is touched; the full arrays are never loaded into RAM. Validates the
    expected dimensionality of each array (chunks >= 1-D, references 2-D,
    reference_lengths 1-D) and confirms the leading-axis row counts agree
    across the triplet.

    Raises:
        MergeError: On any dimensionality or row-count mismatch.
    """
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
    """Fold NumPy's type-promotion rules over a dtype sequence.

    Uses ``np.result_type`` pairwise so the returned dtype can losslessly
    represent every value under NumPy's standard casting rules (e.g. two
    inputs with ``float16`` and ``float32`` promote to ``float32``).
    """
    common = dtypes[0]
    for dtype in dtypes[1:]:
        common = np.result_type(common, dtype)
    return np.dtype(common)


def build_merge_plan(input_directories: Sequence[Path]) -> MergePlan:
    """Inspect every input and compute the merged output layout.

    Reconciles the per-input metadata into global totals: sums row counts,
    checks that every dataset shares the same ``chunks`` tail shape, takes
    the maximum ``references`` width (narrower inputs are later zero-padded),
    and promotes dtypes via ``np.result_type``.

    Raises:
        MergeError: If no inputs are provided, an input path is not a
            directory, or the datasets disagree on ``chunks`` tail shape.
    """
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
    """Return the ``OutputPaths`` triple inside ``output_directory``.

    Uses the same filenames as Bonito (``REQUIRED_FILENAMES``) so the
    merged dataset is a drop-in replacement for a single training split.
    """
    return OutputPaths(
        chunks=output_directory / REQUIRED_FILENAMES[0],
        references=output_directory / REQUIRED_FILENAMES[1],
        reference_lengths=output_directory / REQUIRED_FILENAMES[2],
    )


def ensure_output_paths_are_writable(output_paths: OutputPaths, overwrite: bool) -> None:
    """Guard against clobbering an existing merged dataset unless ``--overwrite``.

    Raises:
        MergeError: If any output file exists and ``overwrite`` is ``False``.
    """
    for output_path in (output_paths.chunks, output_paths.references, output_paths.reference_lengths):
        if output_path.exists() and not overwrite:
            raise MergeError(
                f"Output file already exists: {output_path}. Use --overwrite to replace it."
            )


def create_output_memmaps(plan: MergePlan, output_paths: OutputPaths) -> tuple[np.memmap, np.memmap, np.memmap]:
    """Allocate the three output arrays on disk as writable memmaps.

    Opening with ``mode="w+"`` creates each ``.npy`` file at full final
    size, backed by a sparse file where supported; the OS zeroes pages on
    first touch so no explicit fill is required for ``chunks`` or
    ``reference_lengths``. ``merged_references`` is explicitly zeroed
    afterwards because per-dataset writes only touch the first
    ``reference_width`` columns and the remainder must stay zero to serve
    as the pad region for narrower inputs.
    """
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

    # Explicit zero-fill of references so any trailing columns past a
    # dataset's reference_width remain zero-padded in the merged array.
    merged_references[:] = 0
    return merged_chunks, merged_references, merged_reference_lengths


def estimate_rows_per_block(
    metadata: DatasetMetadata,
    plan: MergePlan,
    target_bytes: int,
) -> int:
    """Size each copy block so the three-array row combined fits in ``target_bytes``.

    Sums the per-row byte cost across all three arrays in the OUTPUT
    dtypes (which may be wider than the inputs after promotion) and floors
    to the target budget. The minimum is clamped to one row so extremely
    large single-row payloads still make progress. The maximum is clamped
    to the dataset's row count so the block loop terminates cleanly on the
    final partial block.
    """
    chunk_elements_per_row = int(math.prod(metadata.chunk_tail_shape)) if metadata.chunk_tail_shape else 1
    chunk_bytes_per_row = chunk_elements_per_row * int(np.dtype(plan.chunks_dtype).itemsize)
    reference_bytes_per_row = metadata.reference_width * int(np.dtype(plan.references_dtype).itemsize)
    reference_length_bytes_per_row = int(np.dtype(plan.reference_lengths_dtype).itemsize)

    bytes_per_row = chunk_bytes_per_row + reference_bytes_per_row + reference_length_bytes_per_row
    # ``max(1, ...)`` on bytes_per_row guards the degenerate zero-cost case.
    return max(1, min(metadata.row_count, target_bytes // max(1, bytes_per_row)))


def stream_merge(
    plan: MergePlan,
    merged_chunks: np.memmap,
    merged_references: np.memmap,
    merged_reference_lengths: np.memmap,
    target_chunk_mb: float,
) -> None:
    """Stream every input triplet into the merged output memmaps block-by-block.

    Datasets are written in the order supplied on the CLI: dataset 0
    occupies rows ``[0, n0)``, dataset 1 rows ``[n0, n0+n1)``, and so on.
    Within each dataset the rows are copied in contiguous blocks sized by
    ``estimate_rows_per_block`` so peak RAM use is bounded by
    ``target_chunk_mb``. References are written into the first
    ``metadata.reference_width`` columns, leaving the remaining columns at
    their zero-initialised state so narrower inputs are effectively
    zero-padded out to the merged width. A final ``flush`` on all three
    memmaps forces dirty pages to disk before the function returns.
    """
    print("[INFO] Streaming arrays into output (mmap -> mmap) ...")
    # MiB -> bytes; ``max(1, ...)`` avoids a zero budget if the user passes
    # a tiny floating-point value.
    target_bytes = max(1, int(target_chunk_mb * 1024 * 1024))
    progress = progress_bar(total=plan.total_rows, desc="Streaming rows", unit="row")

    # Rolling destination cursor: advances exactly by the number of rows
    # copied in each block, independently of the per-dataset cursor below.
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

                # Slice-assignment between two memmaps uses NumPy's
                # implicit dtype cast; the promoted output dtype ensures
                # this stays lossless.
                merged_chunks[destination_row_start:destination_row_end, ...] = source_chunks[
                    source_row_start:source_row_end, ...
                ]
                # Restrict the column slice to the source dataset's own
                # reference_width so the zero pad on the right is
                # preserved for narrower inputs.
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

    # Force dirty mmap pages to disk so a subsequent reopen-for-shuffle
    # observes the full, consistent streamed output.
    merged_chunks.flush()
    merged_references.flush()
    merged_reference_lengths.flush()


def shuffle_output_rows(
    plan: MergePlan,
    output_paths: OutputPaths,
    target_chunk_mb: float,
    seed: int | None,
) -> None:
    """Apply a single row permutation to all three merged arrays on disk.

    The three arrays must stay row-aligned after shuffling, so the same
    permutation is applied in lockstep. Because fancy indexing from a
    memmap back into itself cannot overlap safely, the shuffled output is
    written into ``*.tmp.npy`` memmaps alongside the originals and then
    atomically ``os.replace``d. ``seed`` is forwarded to
    ``np.random.default_rng`` so callers can reproduce a specific shuffle.
    """
    print("[INFO] Shuffling merged rows ...")
    merged_chunks = np.load(output_paths.chunks, allow_pickle=False, mmap_mode="r")
    merged_references = np.load(output_paths.references, allow_pickle=False, mmap_mode="r")
    merged_reference_lengths = np.load(output_paths.reference_lengths, allow_pickle=False, mmap_mode="r")

    # A single permutation over ``total_rows`` guarantees the three
    # arrays remain row-aligned after shuffling.
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

    # Same bytes-per-row accounting as the streaming pass, but uses the
    # merged plan's own dtypes/shapes because no per-dataset metadata is
    # relevant during shuffle.
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
            # Fancy-index with the block's slice of the global permutation.
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

    # ``os.replace`` is atomic on POSIX/Windows, so a crash after a
    # successful rename still leaves a consistent final output.
    os.replace(temporary_paths.chunks, output_paths.chunks)
    os.replace(temporary_paths.references, output_paths.references)
    os.replace(temporary_paths.reference_lengths, output_paths.reference_lengths)


def print_summary(plan: MergePlan, output_paths: OutputPaths, shuffled: bool) -> None:
    """Emit a single ``[DONE]`` block describing shapes, dtypes, and row total.

    Printed to stdout so operators can immediately confirm the merged
    dataset matches their expectation without rerunning ``ReadNumpy.py``.
    """
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
    """CLI entry point: plan, allocate, stream, optional shuffle, and summarise.

    Returns ``0`` on success; ``MergeError``s escape to the top-level
    ``__main__`` block, which prints a concise message and exits with 2.
    """
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
