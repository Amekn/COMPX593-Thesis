#!/usr/bin/env python3
"""Inject the POD5 run-info fields Dorado/Bonito require but simulators omit.

Simulator-generated POD5 files (used in the COMPX593 thesis to build a
controlled deep-mutational-scanning training corpus) ship with sparse
``run_info`` metadata. Dorado refuses to basecall such reads because fields
such as ``flow_cell_product_code``, ``sample_id``, ``run_id`` and the
``context_tags[sample_frequency]`` entry are missing or empty. This script
reads every record, patches those fields in place on a mutable copy, and
writes a new POD5 file; the raw signal payload is not modified.

A single ``run_id`` UUID is generated once per invocation and shared across
every read in the output file, so downstream tooling treats the file as a
single coherent sequencing run.
"""

from __future__ import annotations

import argparse
import uuid
from dataclasses import replace
from pathlib import Path
from typing import Sequence

import pod5

# MinION R10.4.1 flow cell code. Matches the chemistry Dorado's built-in
# models are trained against, so basecalling runs without a model override.
DEFAULT_FLOWCELL_CODE = "FLO-MIN114"
# Sentinel sample name flagging records that originated from simulation
# rather than a real sequencing experiment.
DEFAULT_SAMPLE_ID = "SIMULATED"
# 5 kHz is the current-generation MinION/PromethION sampling rate used
# throughout the thesis; must match the model's expected cadence.
DEFAULT_SAMPLE_FREQUENCY = 5000


def build_argument_parser() -> argparse.ArgumentParser:
    """Create the CLI parser for the metadata injector.

    All metadata values default to the constants above so the common case
    (simulated 5 kHz R10.4.1 reads) needs no flags. Overrides are exposed
    for flow cell code, sample id, and sampling frequency.
    """
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
    """Validate the injection request before any output file is touched.

    Raises:
        FileNotFoundError: If the input POD5 does not exist.
        ValueError: If extensions are wrong or ``sample_frequency`` is not
            a positive integer.
        FileExistsError: If the output exists and ``overwrite`` was not set.
    """
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
    """Rewrite ``read.run_info`` with Dorado-compatible metadata in place.

    ``pod5.RunInfo`` is a frozen dataclass, so the function builds a new
    instance via ``dataclasses.replace`` and assigns it back. The caller is
    responsible for supplying a single ``run_id`` that stays constant across
    every read of the output file, otherwise Dorado will treat the reads as
    originating from unrelated runs.

    Args:
        read: Mutable read obtained from ``ReadRecord.to_read()``.
        flowcell_code: Value stored in ``flow_cell_product_code``.
        sample_id: Value stored in ``sample_id``.
        sample_frequency: Sampling rate in Hz; stored as a string in
            ``context_tags['sample_frequency']`` because POD5 persists
            context tags as ``str -> str`` dicts.
        run_id: Shared run UUID string.
    """
    # ``tracking_id`` and ``context_tags`` may be None on simulated reads;
    # copy into a new dict so the frozen RunInfo's backing map is untouched.
    tracking_id = dict(getattr(read.run_info, "tracking_id", {}) or {})
    # ``exp_start_time`` is required by some Dorado versions; populate with
    # an epoch-equivalent sentinel only when absent to preserve real values.
    tracking_id.setdefault("exp_start_time", "0.0")
    tracking_id["run_id"] = run_id

    context_tags = dict(getattr(read.run_info, "context_tags", {}) or {})
    # Context tags are persisted as strings; cast explicitly to avoid a
    # schema mismatch when the POD5 writer serialises the map.
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
    """Stream every read of the input, patch metadata, and write a new POD5.

    The function generates exactly one run UUID per invocation and shares
    it across every read so the output POD5 reads back as a single run.

    Args:
        input_path: Source POD5 file.
        output_path: Destination POD5 file (must be writable; the caller is
            responsible for honouring ``--overwrite`` semantics).
        flowcell_code: Flow cell product code written into run_info.
        sample_id: Sample identifier written into run_info.
        sample_frequency: Sampling frequency in Hz written into
            ``run_info.context_tags``.

    Returns:
        The number of reads processed and written.
    """
    shared_run_id = str(uuid.uuid4())
    processed_read_count = 0

    with pod5.Reader(str(input_path)) as reader, pod5.Writer(str(output_path)) as writer:
        for read_record in reader:
            # ``to_read()`` yields a detached, mutable Read. Mutating the
            # original ReadRecord from the reader is not supported.
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
    """CLI entry point: validate arguments, run injection, print count."""
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
