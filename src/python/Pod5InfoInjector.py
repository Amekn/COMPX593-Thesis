#!/usr/bin/env python3
"""Populate missing POD5 run metadata needed by downstream basecallers."""

from __future__ import annotations

import argparse
import uuid
from dataclasses import replace
from pathlib import Path
from typing import Sequence

import pod5

DEFAULT_FLOWCELL_CODE = "FLO-MIN114"
DEFAULT_SAMPLE_ID = "SIMULATED"
DEFAULT_SAMPLE_FREQUENCY = 5000


def build_argument_parser() -> argparse.ArgumentParser:
    """Create the command-line interface for metadata injection."""
    parser = argparse.ArgumentParser(
        description=(
            "Read a POD5 file, populate required run_info fields for Dorado compatibility, "
            "and write the updated reads to a new POD5 file."
        )
    )
    parser.add_argument("input_pod5", type=Path, help="Source POD5 file.")
    parser.add_argument("output_pod5", type=Path, help="Destination POD5 file.")
    parser.add_argument(
        "--flowcell-code",
        default=DEFAULT_FLOWCELL_CODE,
        help="Flow cell product code to store in run_info.",
    )
    parser.add_argument(
        "--sample-id",
        default=DEFAULT_SAMPLE_ID,
        help="Sample identifier to store in run_info.",
    )
    parser.add_argument(
        "--sample-frequency",
        type=int,
        default=DEFAULT_SAMPLE_FREQUENCY,
        help="Sampling frequency, in Hz, to store in run_info.context_tags.",
    )
    parser.add_argument(
        "--overwrite",
        action="store_true",
        help="Allow replacing an existing output file.",
    )
    return parser


def validate_arguments(input_path: Path, output_path: Path, overwrite: bool, sample_frequency: int) -> None:
    """Validate the metadata injection request before any output is written."""
    if not input_path.exists():
        raise FileNotFoundError(f"Input POD5 file does not exist: {input_path}")
    if input_path.suffix.lower() != ".pod5":
        raise ValueError(f"Input path must use the .pod5 extension: {input_path}")
    if output_path.suffix.lower() != ".pod5":
        raise ValueError(f"Output path must use the .pod5 extension: {output_path}")
    if output_path.exists() and not overwrite:
        raise FileExistsError(f"Output file already exists: {output_path}")
    if sample_frequency <= 0:
        raise ValueError("sample-frequency must be a positive integer")


def populate_run_info(
    read: pod5.Read,
    flowcell_code: str,
    sample_id: str,
    sample_frequency: int,
    run_id: str,
) -> None:
    """Update the mutable POD5 read with run metadata required by downstream tools."""
    tracking_id = dict(getattr(read.run_info, "tracking_id", {}) or {})
    tracking_id.setdefault("exp_start_time", "0.0")
    tracking_id["run_id"] = run_id

    context_tags = dict(getattr(read.run_info, "context_tags", {}) or {})
    context_tags["sample_frequency"] = str(sample_frequency)

    read.run_info = replace(
        read.run_info,
        flow_cell_product_code=flowcell_code,
        sample_id=sample_id,
        tracking_id=tracking_id,
        context_tags=context_tags,
    )


def inject_run_info(
    input_path: Path,
    output_path: Path,
    flowcell_code: str,
    sample_id: str,
    sample_frequency: int,
) -> int:
    """Write a new POD5 file with populated run metadata and return the read count."""
    shared_run_id = str(uuid.uuid4())
    processed_read_count = 0

    with pod5.Reader(str(input_path)) as reader, pod5.Writer(str(output_path)) as writer:
        for read_record in reader:
            mutable_read = read_record.to_read()
            populate_run_info(
                read=mutable_read,
                flowcell_code=flowcell_code,
                sample_id=sample_id,
                sample_frequency=sample_frequency,
                run_id=shared_run_id,
            )
            writer.add_read(mutable_read)
            processed_read_count += 1

    return processed_read_count


def main(argv: Sequence[str] | None = None) -> int:
    """Run the command-line entry point."""
    arguments = build_argument_parser().parse_args(argv)

    input_path = arguments.input_pod5.resolve()
    output_path = arguments.output_pod5.resolve()
    validate_arguments(
        input_path=input_path,
        output_path=output_path,
        overwrite=arguments.overwrite,
        sample_frequency=arguments.sample_frequency,
    )
    output_path.parent.mkdir(parents=True, exist_ok=True)

    processed_read_count = inject_run_info(
        input_path=input_path,
        output_path=output_path,
        flowcell_code=arguments.flowcell_code,
        sample_id=arguments.sample_id,
        sample_frequency=arguments.sample_frequency,
    )
    print(f"Updated {processed_read_count} reads and wrote {output_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
