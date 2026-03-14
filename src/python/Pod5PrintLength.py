#!/usr/bin/env python3
"""Print the signal length for every read in a POD5 file."""

from __future__ import annotations

import argparse
from pathlib import Path
from typing import Sequence

import pod5


def build_argument_parser() -> argparse.ArgumentParser:
    """Create the command-line interface for reporting POD5 read lengths."""
    parser = argparse.ArgumentParser(
        description="Print '<read_id>\\t<signal_length>' for every read in a POD5 file."
    )
    parser.add_argument("pod5_path", type=Path, help="Input POD5 file.")
    return parser


def main(argv: Sequence[str] | None = None) -> int:
    """Run the command-line entry point."""
    arguments = build_argument_parser().parse_args(argv)
    pod5_path = arguments.pod5_path.resolve()

    if not pod5_path.exists():
        raise FileNotFoundError(f"Input POD5 file does not exist: {pod5_path}")
    if pod5_path.suffix.lower() != ".pod5":
        raise ValueError(f"Input path must use the .pod5 extension: {pod5_path}")

    with pod5.Reader(str(pod5_path)) as reader:
        for read_record in reader.reads():
            signal_length = int(getattr(read_record, "num_samples", len(read_record.signal)))
            print(f"{read_record.read_id}\t{signal_length}")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
