#!/usr/bin/env python3
"""Render a single ONT raw-signal trace from a POD5 file as a publication SVG.

Produces the raw current-vs-time figure used in the COMPX593 thesis to
illustrate a single nanopore read. Given a POD5 path and a UUID read_id,
the script selects the matching read via the official pod5 reader API,
converts the signal to picoamps (preferring the cached ``signal_pa`` view
and falling back to on-the-fly calibration of the raw ADC array), builds a
minimalist serif-styled figure, optionally downsamples to keep the SVG
lightweight, and writes the result as SVG. Exit codes are used to
distinguish invalid UUIDs (2), missing read_ids (3), missing POD5 files (1),
and unexpected errors (99).

Usage:
    python PlotSignal.py <file.pod5> <read_id_uuid> [--output figure.svg]

Dependencies:
    pip install pod5 numpy matplotlib
"""

import argparse
from pathlib import Path
import sys
import uuid
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.figure import Figure
from matplotlib.ticker import AutoMinorLocator

try:
    import pod5 as p5
except Exception as e:
    print("Error: the 'pod5' package is required. Install with 'pip install pod5'.", file=sys.stderr)
    raise


def configure_matplotlib() -> None:
    """Apply the thesis-wide matplotlib rcParams for figure consistency.

    Settings target a serif-style journal aesthetic: DejaVu Serif body text,
    thin axes, no top/right spines, no grid by default, and ``svg.fonttype``
    set to ``"none"`` so text remains editable/selectable in the vector
    output rather than being rasterised into glyph paths.
    """
    plt.rcParams.update(
        {
            "font.family": "DejaVu Serif",
            "font.size": 10,
            "axes.titlesize": 11,
            "axes.labelsize": 10,
            "axes.linewidth": 0.8,
            "xtick.labelsize": 9,
            "ytick.labelsize": 9,
            "xtick.direction": "out",
            "ytick.direction": "out",
            "xtick.major.width": 0.8,
            "ytick.major.width": 0.8,
            "xtick.minor.width": 0.6,
            "ytick.minor.width": 0.6,
            "xtick.major.size": 4,
            "ytick.major.size": 4,
            "xtick.minor.size": 2.5,
            "ytick.minor.size": 2.5,
            "figure.facecolor": "white",
            "axes.facecolor": "white",
            "axes.edgecolor": "#222222",
            "axes.grid": False,
            "savefig.facecolor": "white",
            "savefig.edgecolor": "white",
            "savefig.bbox": "tight",
            "savefig.pad_inches": 0.03,
            "svg.fonttype": "none",
        }
    )


def maybe_downsample(time_s: np.ndarray, signal_pa: np.ndarray, max_points: int) -> tuple[np.ndarray, np.ndarray]:
    """Uniformly subsample ``(time, signal)`` to at most ``max_points`` samples.

    Long nanopore reads at 5 kHz can exceed millions of samples; a plain SVG
    plot of that many path vertices bloats the file and is slower to render
    than it is to acquire. A uniform stride via ``np.linspace`` preserves
    the overall shape well enough for a whole-read overview figure. Pass
    ``max_points <= 0`` to disable downsampling and plot every sample.

    Returns:
        The possibly-subsampled ``(time_s, signal_pa)`` pair. No copy is
        made when the input already fits.
    """
    if max_points <= 0 or signal_pa.size <= max_points:
        return time_s, signal_pa

    # int64 indexing avoids overflow for very long reads on 32-bit builds.
    idx = np.linspace(0, signal_pa.size - 1, max_points, dtype=np.int64)
    return time_s[idx], signal_pa[idx]


def build_output_path(pod5_path: str, read_id: str, output: str | None) -> Path:
    """Resolve the SVG destination, defaulting beside the source POD5.

    If ``output`` is supplied, it is used verbatim; otherwise the file is
    placed next to the POD5 as ``<pod5_stem>_<read_id>.svg``. The suffix is
    normalised to ``.svg`` because the script only writes SVG.
    """
    if output:
        out_path = Path(output)
    else:
        out_path = Path(pod5_path).with_name(f"{Path(pod5_path).stem}_{read_id}.svg")

    if out_path.suffix.lower() != ".svg":
        out_path = out_path.with_suffix(".svg")

    return out_path


def make_figure(
    time_s: np.ndarray,
    signal_pa: np.ndarray,
    read_id: str,
    sample_rate: float,
    n_samples: int,
    channel: object,
) -> Figure:
    """Build the single-axes current-vs-time figure for the thesis.

    The figure is laid out at a 7.2x2.8-inch landscape aspect ratio (close
    to the typical single-column body-text width used in the thesis). The
    trace is drawn as a thin navy line; axis metadata (read id, sample
    rate, sample count, optional channel) is embedded in the title for
    self-contained figure captions. Top and right spines are hidden and a
    5% vertical pad is added so the trace peaks aren't clipped.
    """
    fig, ax = plt.subplots(figsize=(7.2, 2.8), constrained_layout=True)
    ax.plot(
        time_s,
        signal_pa,
        color="#0f4c81",
        linewidth=0.75,
        alpha=0.95,
        solid_capstyle="round",
        solid_joinstyle="round",
        rasterized=False,
    )

    ax.set_xlabel("Time (s)")
    ax.set_ylabel("Current (pA)")

    meta = [f"id = {read_id}", f"s = {int(sample_rate)} Hz", f"n = {n_samples:,}"]
    if channel is not None:
        meta.insert(1, f"c = {channel}")
    ax.set_title(" | ".join(meta), loc="left", pad=8)

    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.xaxis.set_minor_locator(AutoMinorLocator(2))
    ax.yaxis.set_minor_locator(AutoMinorLocator(2))
    ax.grid(axis="y", color="#d9d9d9", linewidth=0.6, alpha=0.7)
    # Collapse the x-axis padding so the trace starts flush at t=0.
    ax.margins(x=0)

    # Explicit 5% vertical pad prevents the plotted line from touching the
    # axes frame; the guard against a zero-range signal avoids
    # ``set_ylim(0, 0)`` which matplotlib rejects.
    y_min = float(np.min(signal_pa))
    y_max = float(np.max(signal_pa))
    if y_max > y_min:
        pad = 0.05 * (y_max - y_min)
        ax.set_ylim(y_min - pad, y_max + pad)

    return fig

def main() -> int:
    """CLI entry point: parse args, locate the read, render, and save SVG.

    Returns:
        ``0`` on success; ``1`` if the POD5 path cannot be opened; ``2`` if
        the supplied read_id is not a valid UUID; ``3`` if the UUID is not
        present in the POD5; ``99`` for any other unhandled exception.
    """
    configure_matplotlib()

    parser = argparse.ArgumentParser(
        description="Plot signal vs time for a specific read_id from a POD5 file and save as SVG."
    )
    parser.add_argument("pod5_path", help="Path to the .pod5 file")
    parser.add_argument("read_id", help="UUID read_id to plot (e.g. 0000173c-bf67-...)")
    parser.add_argument(
        "-o",
        "--output",
        help="Output SVG path. Defaults to <pod5_stem>_<read_id>.svg next to the POD5 file.",
    )
    parser.add_argument(
        "--max-points",
        type=int,
        default=20000,
        help="Maximum plotted points for SVG output. Use 0 to disable downsampling. Default: 20000.",
    )
    parser.add_argument(
        "--show",
        action="store_true",
        help="Also display the figure interactively after saving.",
    )
    args = parser.parse_args()

    # Validate and normalize the UUID format early to avoid surprises.
    try:
        read_id = str(uuid.UUID(args.read_id))
    except ValueError:
        print(f"Error: '{args.read_id}' is not a valid UUID read_id.", file=sys.stderr)
        return 2

    try:
        with p5.Reader(args.pod5_path) as reader:  # official Reader API
            # Reader.reads() accepts an iterable of read_id strings (UUIDs).
            # next(...) raises StopIteration if not found.
            try:
                read = next(reader.reads([read_id]))
            except StopIteration:
                print(
                    f"Error: read_id {read_id} not found in {args.pod5_path}.\n"
                    "Tip: list IDs quickly with: pod5 view /path/to/file.pod5 --include \"read_id\" | head",
                    file=sys.stderr,
                )
                return 3

            # Sample rate (Hz) for the run; used to create the time axis.
            sample_rate = float(read.run_info.sample_rate)  # Hz

            # Prefer calibrated picoamp signal if available; otherwise calibrate raw ADC.
            # Both interfaces are part of the official API.
            if hasattr(read, "signal_pa"):
                signal_pa = read.signal_pa
            else:
                # Fallback for older versions: calibrate the int16 ADC array.
                signal_pa = read.calibrate_signal_array(read.signal)

            n = int(read.num_samples)
            t = np.arange(n, dtype=np.float64) / sample_rate
            plot_t, plot_signal = maybe_downsample(t, np.asarray(signal_pa), args.max_points)

            ch = getattr(getattr(read, "pore", None), "channel", None)
            out_path = build_output_path(args.pod5_path, read_id, args.output)
            fig = make_figure(plot_t, plot_signal, read_id, sample_rate, n, ch)
            out_path.parent.mkdir(parents=True, exist_ok=True)
            fig.savefig(out_path, format="svg")
            print(f"Saved SVG plot to: {out_path}")

            if args.show:
                plt.show()
            else:
                plt.close(fig)

    except FileNotFoundError:
        print(f"Error: file not found: {args.pod5_path}", file=sys.stderr)
        return 1
    except Exception as e:
        print(f"Unhandled error: {e}", file=sys.stderr)
        return 99

    return 0

if __name__ == "__main__":
    raise SystemExit(main())
