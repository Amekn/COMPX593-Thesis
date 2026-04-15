#!/usr/bin/env python3
"""Print the raw-signal length (in samples) of every read in a POD5 file.

Diagnostic tool used during POD5 preprocessing for the COMPX593 thesis to
audit per-read sample counts before feeding data to Dorado/Bonito. One
tab-separated line per read (``<read_id>\\t<signal_length>``) is emitted to
stdout so the output composes with standard shell pipelines (``wc -l``,
``awk``, ``sort``). No files are modified.
"""

from __future__ import annotations

import argparse
from pathlib import Path
from typing import Sequence

import pod5


def build_argument_parser() -> argparse.ArgumentParser:
    """Create the CLI parser accepting a single POD5 input path."""
    parser = argparse.ArgumentParser(
        description="Print '<read_id>\\t<signal_length>' for every read in a POD5 file."
    )
    parser.add_argument("pod5_path", type=Path, help="Input POD5 file.")
    return parser


def main(argv: Sequence[str] | None = None) -> int:
    """Stream every read in the POD5 file and print its signal length.

    Args:
        argv: Optional argument vector for testing; defaults to ``sys.argv[1:]``.

    Returns:
        Process exit status (always ``0`` on success).

    Raises:
        FileNotFoundError: If the POD5 path does not exist.
        ValueError: If the path does not carry the ``.pod5`` extension.
    """
    arguments = build_argument_parser().parse_args(argv)
    pod5_path = arguments.pod5_path.resolve()

    if not pod5_path.exists():
        raise FileNotFoundError(f"Input POD5 file does not exist: {pod5_path}")
    if pod5_path.suffix.lower() != ".pod5":
        raise ValueError(f"Input path must use the .pod5 extension: {pod5_path}")

    with pod5.Reader(str(pod5_path)) as reader:
        for read_record in reader.reads():
            # Prefer ``num_samples`` (O(1) header read) and fall back to the
            # materialised signal length for older pod5 API versions where
            # ``num_samples`` is not exposed on the read record.
            signal_length = int(getattr(read_record, "num_samples", len(read_record.signal)))
            print(f"{read_record.read_id}\t{signal_length}")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
