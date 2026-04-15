#!/usr/bin/env python3
"""Report the shape of a NumPy array stored in a ``.npy`` file.

Diagnostic helper used throughout the COMPX593 thesis pipeline to verify
Bonito training-array dimensions (e.g. leading-axis row counts of
``chunks.npy`` / ``references.npy`` / ``reference_lengths.npy``) without
deserialising the full array into memory beyond what ``np.load`` requires.
Prints the shape tuple to stdout; returns exit status only.
"""

from __future__ import annotations

import argparse
from pathlib import Path
from typing import Sequence

import numpy as np


def build_argument_parser() -> argparse.ArgumentParser:
    """Create the CLI parser that accepts a single ``.npy`` path argument."""
    parser = argparse.ArgumentParser(
        description="Load a .npy array file and print its shape."
    )
    parser.add_argument("input_numpy_file", type=Path, help="Input .npy file.")
    return parser


def main(argv: Sequence[str] | None = None) -> int:
    """Load the target ``.npy`` file and print its shape tuple.

    Args:
        argv: Optional argument vector for testing; defaults to ``sys.argv[1:]``.

    Returns:
        Process exit status (always ``0`` on success).

    Raises:
        FileNotFoundError: If the supplied path does not exist.
        ValueError: If the path does not carry the ``.npy`` extension.
    """
    arguments = build_argument_parser().parse_args(argv)
    numpy_path = arguments.input_numpy_file.resolve()

    if not numpy_path.exists():
        raise FileNotFoundError(f"Input array file does not exist: {numpy_path}")
    if numpy_path.suffix.lower() != ".npy":
        raise ValueError(f"Input path must use the .npy extension: {numpy_path}")

    # allow_pickle=False guards against malicious object arrays; Bonito training
    # arrays are plain numeric tensors, so this is safe and faster.
    array = np.load(numpy_path, allow_pickle=False)
    print(array.shape)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
