#!/usr/bin/env python3
"""Report the shape of a NumPy array stored in a .npy file."""

from __future__ import annotations

import argparse
from pathlib import Path
from typing import Sequence

import numpy as np


def build_argument_parser() -> argparse.ArgumentParser:
    """Create the command-line interface for array inspection."""
    parser = argparse.ArgumentParser(
        description="Load a .npy array file and print its shape."
    )
    parser.add_argument("input_numpy_file", type=Path, help="Input .npy file.")
    return parser


def main(argv: Sequence[str] | None = None) -> int:
    """Run the command-line entry point."""
    arguments = build_argument_parser().parse_args(argv)
    numpy_path = arguments.input_numpy_file.resolve()

    if not numpy_path.exists():
        raise FileNotFoundError(f"Input array file does not exist: {numpy_path}")
    if numpy_path.suffix.lower() != ".npy":
        raise ValueError(f"Input path must use the .npy extension: {numpy_path}")

    array = np.load(numpy_path, allow_pickle=False)
    print(array.shape)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
