#!/usr/bin/env python3
"""
Merge Bonito-style training arrays from multiple folders, **without** loading
all arrays into RAM at once, with progress bars and ETA.

Strategy (RAM-safe):
  1) First pass ("Inspecting headers"): open each input .npy with
     mmap_mode='r' to read only shapes/dtypes and count rows; track max
     width for references.
  2) Create output .npy as memory-mapped arrays using
     numpy.lib.format.open_memmap.
  3) Second pass ("Streaming rows"): copy in **row blocks** to bound RAM and
     show smooth progress with ETA. For references.npy, copy only actual width
     per folder into the left columns; output is zero-initialized for padding.

Each input folder must contain:
  - chunks.npy
  - references.npy
  - reference_lengths.npy

Behavior:
  - chunks.npy: concatenated along axis=0.
  - references.npy: right-zero-padded to the maximum width across inputs,
                    then concatenated along axis=0.
  - reference_lengths.npy: concatenated along axis=0.

Example:
    python merge_numpy.py -o merged_dir input1 input2 input3 --overwrite \
        --target-chunk-mb 64
"""

from __future__ import annotations
import argparse
from pathlib import Path
import sys
import os
from typing import List, Tuple
import math
import numpy as np
from numpy.lib.format import open_memmap

# Optional tqdm with a safe fallback (no-op with same API used here)
try:
    from tqdm import tqdm as _tqdm
except Exception:  # pragma: no cover
    class _TqdmNoOp:
        def __init__(self, total=None, desc=None, unit=None, **kwargs):
            self.total = total
        def update(self, n=1):
            pass
        def set_postfix(self, **kwargs):
            pass
        def close(self):
            pass
    def _tqdm(*args, **kwargs):  # type: ignore
        return _TqdmNoOp(*args, **kwargs)

REQ_FILES = ("chunks.npy", "references.npy", "reference_lengths.npy")

class MergeError(Exception):
    pass

# ──────────────────────────────────────────────────────────────────────────────
# Helpers
# ──────────────────────────────────────────────────────────────────────────────

def check_required_files(folder: Path) -> Tuple[Path, Path, Path]:
    chunks = folder / REQ_FILES[0]
    refs = folder / REQ_FILES[1]
    lengths = folder / REQ_FILES[2]

    missing = [str(p.name) for p in (chunks, refs, lengths) if not p.exists()]
    if missing:
        raise MergeError(f"Missing files in {folder}: {', '.join(missing)}")
    return chunks, refs, lengths


def inspect_triplet(folder: Path):
    """Open arrays in mmap read-only mode and return lightweight metadata.

    Returns:
        (n_rows, chunks_tail_shape, refs_width, dtypes_tuple)
    where dtypes_tuple is (chunks_dtype, refs_dtype, reflens_dtype).
    """
    chunks_p, refs_p, lengths_p = check_required_files(folder)

    chunks = np.load(chunks_p, allow_pickle=False, mmap_mode='r')
    refs = np.load(refs_p, allow_pickle=False, mmap_mode='r')
    reflens = np.load(lengths_p, allow_pickle=False, mmap_mode='r')

    if chunks.ndim < 1:
        raise MergeError(f"chunks.npy must be at least 1-D; got {chunks.shape}")
    if refs.ndim != 2:
        raise MergeError(f"references.npy must be 2-D; got {refs.shape}")
    if reflens.ndim != 1:
        raise MergeError(f"reference_lengths.npy must be 1-D; got {reflens.shape}")

    n_chunks = chunks.shape[0]
    n_refs = refs.shape[0]
    n_lens = reflens.shape[0]
    if not (n_chunks == n_refs == n_lens):
        raise MergeError(
            f"Row count mismatch in {folder} (chunks: {n_chunks}, references: {n_refs}, reference_lengths: {n_lens})"
        )

    chunks_tail_shape = chunks.shape[1:]  # may be (T,) or (T, F, ...)
    refs_width = refs.shape[1]
    dtypes = (chunks.dtype, refs.dtype, reflens.dtype)

    return n_chunks, chunks_tail_shape, refs_width, dtypes


def common_dtype(dtypes: List[np.dtype]) -> np.dtype:
    """Compute a safe common dtype via np.result_type."""
    dt = dtypes[0]
    for d in dtypes[1:]:
        dt = np.result_type(dt, d)
    return np.dtype(dt)

# ──────────────────────────────────────────────────────────────────────────────
# CLI
# ──────────────────────────────────────────────────────────────────────────────

def parse_args():
    parser = argparse.ArgumentParser(
        description=(
            "Merge Bonito training arrays (chunks/references/reference_lengths) "
            "from multiple folders without high RAM usage, with progress bars."
        )
    )
    parser.add_argument(
        "input_dirs",
        nargs="+",
        help="Input folders (order is preserved). Each must contain the three .npy files.",
    )
    parser.add_argument(
        "-o",
        "--output",
        required=True,
        help="Output folder to write merged arrays into.",
    )
    parser.add_argument(
        "--overwrite",
        action="store_true",
        help="Allow overwriting existing output files.",
    )
    parser.add_argument(
        "--target-chunk-mb",
        type=float,
        default=64.0,
        help=(
            "Approximate memory budget per streaming block (in MiB). "
            "Larger values can be faster but use more RAM."
        ),
    )
    parser.add_argument(
        "--shuffle",
        action="store_true",
        help="Shuffle rows randomly after merge, keeping correspondence among arrays.",
    )
    parser.add_argument(
        "--seed",
        type=int,
        default=None,
        help="Random seed for shuffling (default: system entropy).",
    )
    return parser.parse_args()

# ──────────────────────────────────────────────────────────────────────────────
# Main
# ──────────────────────────────────────────────────────────────────────────────

def main() -> int:
    args = parse_args()

    input_paths = [Path(p).resolve() for p in args.input_dirs]
    out_dir = Path(args.output).resolve()
    out_dir.mkdir(parents=True, exist_ok=True)

    print(f"[INFO] Loading headers from {len(input_paths)} input folder(s) (mmap, no bulk RAM).")

    # First pass: collect metadata only
    total_rows = 0
    chunk_tail_shapes = set()
    refs_widths: List[int] = []

    chunk_dtypes: List[np.dtype] = []
    refs_dtypes: List[np.dtype] = []
    reflens_dtypes: List[np.dtype] = []

    per_folder_rows: List[int] = []
    per_folder_refs_width: List[int] = []

    bar1 = _tqdm(total=len(input_paths), desc="Inspecting headers", unit="dir")
    for folder in input_paths:
        if not folder.is_dir():
            raise MergeError(f"Not a directory: {folder}")
        n_rows, chunks_tail_shape, refs_width, dtypes = inspect_triplet(folder)

        total_rows += n_rows
        chunk_tail_shapes.add(chunks_tail_shape)
        refs_widths.append(refs_width)

        chunk_dtypes.append(dtypes[0])
        refs_dtypes.append(dtypes[1])
        reflens_dtypes.append(dtypes[2])

        per_folder_rows.append(n_rows)
        per_folder_refs_width.append(refs_width)
        bar1.update(1)
    bar1.close()

    if len(chunk_tail_shapes) != 1:
        detail = ", ".join(sorted(map(str, chunk_tail_shapes)))
        raise MergeError(f"Inconsistent chunks tail shapes across inputs: {detail}")

    max_ref_w = max(refs_widths) if refs_widths else 0

    chunks_tail_shape = next(iter(chunk_tail_shapes))
    chunks_common_dtype = common_dtype(chunk_dtypes)
    refs_common_dtype = common_dtype(refs_dtypes)
    reflens_common_dtype = common_dtype(reflens_dtypes)

    # Prepare output paths
    out_chunks = out_dir / "chunks.npy"
    out_refs = out_dir / "references.npy"
    out_reflens = out_dir / "reference_lengths.npy"

    # Overwrite checks
    for p in (out_chunks, out_refs, out_reflens):
        if p.exists() and not args.overwrite:
            raise MergeError(f"Output file exists: {p}. Use --overwrite to replace.")

    # Create memmapped outputs
    print("[INFO] Creating output memmaps ...")
    chunks_shape = (total_rows,) + chunks_tail_shape
    refs_shape = (total_rows, max_ref_w)
    reflens_shape = (total_rows,)

    out_chunks_mm = open_memmap(
        out_chunks, mode='w+', dtype=chunks_common_dtype, shape=chunks_shape
    )
    out_refs_mm = open_memmap(
        out_refs, mode='w+', dtype=refs_common_dtype, shape=refs_shape
    )
    out_reflens_mm = open_memmap(
        out_reflens, mode='w+', dtype=reflens_common_dtype, shape=reflens_shape
    )

    # Initialize references with zeros for padding
    out_refs_mm[:] = 0

    # Second pass: stream-copy into slices with progress/ETA
    print("[INFO] Streaming arrays into output (mmap → mmap) ...")

    # Compute per-row byte size to choose a block size near target MiB
    # chunks bytes per row
    chunks_row_elems = int(math.prod(chunks_tail_shape)) if len(chunks_tail_shape) else 1
    chunks_row_bytes = chunks_row_elems * int(np.dtype(chunks_common_dtype).itemsize)
    # refs bytes per row depends on *actual* width per folder; we compute per folder
    refs_itemsize = int(np.dtype(refs_common_dtype).itemsize)
    # reflens is 1 per row
    reflens_row_bytes = int(np.dtype(reflens_common_dtype).itemsize)

    target_bytes = max(1, int(args.target_chunk_mb * 1024 * 1024))

    bar2 = _tqdm(total=total_rows, desc="Streaming rows", unit="row")
    write_pos = 0
    for folder_idx, folder in enumerate(input_paths):
        n_rows = per_folder_rows[folder_idx]
        ref_w = per_folder_refs_width[folder_idx]
        r_end = write_pos + n_rows

        # Estimate bytes per row for this folder
        per_row_bytes = chunks_row_bytes + ref_w * refs_itemsize + reflens_row_bytes
        block_rows = max(1, min(n_rows, target_bytes // max(1, per_row_bytes)))

        # Load inputs as memmaps
        chunks_p, refs_p, lengths_p = check_required_files(folder)
        chunks = np.load(chunks_p, allow_pickle=False, mmap_mode='r')
        refs = np.load(refs_p, allow_pickle=False, mmap_mode='r')
        reflens = np.load(lengths_p, allow_pickle=False, mmap_mode='r')

        # Block-wise copy with progress updates
        start = 0
        while start < n_rows:
            end = min(n_rows, start + block_rows)
            rs, re = write_pos + start, write_pos + end

            # chunks
            out_chunks_mm[rs:re, ...] = chunks[start:end, ...]
            # references (only copy actual width)
            if ref_w > 0:
                out_refs_mm[rs:re, :ref_w] = refs[start:end, :ref_w]
            # reference lengths
            out_reflens_mm[rs:re] = reflens[start:end]

            bar2.update(end - start)
            start = end

        write_pos = r_end
    bar2.close()

    # Ensure data is flushed to disk
    out_chunks_mm.flush()
    out_refs_mm.flush()
    out_reflens_mm.flush()

    # ─────────────────────────────────────────────────────────────────────────
    # Optional shuffling to randomize row order while preserving alignment
    # ─────────────────────────────────────────────────────────────────────────
    if args.shuffle:
        print("[INFO] Shuffling merged rows ...")
        rng = np.random.default_rng(args.seed)
        perm = rng.permutation(total_rows)

        # Estimate per‑row size using the maximum reference width
        per_row_bytes_max = chunks_row_bytes + max_ref_w * refs_itemsize + reflens_row_bytes
        shuffle_block_rows = max(1, int(target_bytes // max(1, per_row_bytes_max)))

        # Temporary files for shuffled output
        tmp_chunks = out_dir / "chunks.tmp.npy"
        tmp_refs = out_dir / "references.tmp.npy"
        tmp_reflens = out_dir / "reference_lengths.tmp.npy"

        tmp_chunks_mm = open_memmap(tmp_chunks, mode='w+', dtype=chunks_common_dtype, shape=chunks_shape)
        tmp_refs_mm = open_memmap(tmp_refs, mode='w+', dtype=refs_common_dtype, shape=refs_shape)
        tmp_reflens_mm = open_memmap(tmp_reflens, mode='w+', dtype=reflens_common_dtype, shape=reflens_shape)

        bar3 = _tqdm(total=total_rows, desc="Shuffling", unit="row")
        for start in range(0, total_rows, shuffle_block_rows):
            end = min(total_rows, start + shuffle_block_rows)
            idx = perm[start:end]
            tmp_chunks_mm[start:end, ...] = out_chunks_mm[idx, ...]
            tmp_refs_mm[start:end, ...] = out_refs_mm[idx, ...]
            tmp_reflens_mm[start:end] = out_reflens_mm[idx]
            bar3.update(end - start)
        bar3.close()

        tmp_chunks_mm.flush()
        tmp_refs_mm.flush()
        tmp_reflens_mm.flush()

        # Atomically replace originals with shuffled versions
        os.replace(tmp_chunks, out_chunks)
        os.replace(tmp_refs, out_refs)
        os.replace(tmp_reflens, out_reflens)

    # Final report
    print(f"[DONE] Wrote merged arrays to {out_dir}:")
    if args.shuffle:
        print("       (rows were shuffled)")
    print(f"       {out_chunks.name}:   shape {chunks_shape}, dtype {out_chunks_mm.dtype}")
    print(f"       {out_refs.name}:     shape {refs_shape}, dtype {out_refs_mm.dtype} (max width {max_ref_w})")
    print(f"       {out_reflens.name}:  shape {reflens_shape}, dtype {out_reflens_mm.dtype}")
    print(f"[DONE] Total records: {total_rows}")
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except MergeError as e:
        print(f"[ERROR] {e}", file=sys.stderr)
        raise SystemExit(2)