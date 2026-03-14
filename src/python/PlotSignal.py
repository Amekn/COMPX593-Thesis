#!/usr/bin/env python3
"""
Plot the ONT raw signal for a single read from a .pod5 file and save a
publication-grade SVG figure.

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
    """Use a clean, journal-friendly plotting style."""
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
    """Limit plotted points so long reads remain sharp and lightweight in SVG."""
    if max_points <= 0 or signal_pa.size <= max_points:
        return time_s, signal_pa

    idx = np.linspace(0, signal_pa.size - 1, max_points, dtype=np.int64)
    return time_s[idx], signal_pa[idx]


def build_output_path(pod5_path: str, read_id: str, output: str | None) -> Path:
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
    ax.margins(x=0)

    y_min = float(np.min(signal_pa))
    y_max = float(np.max(signal_pa))
    if y_max > y_min:
        pad = 0.05 * (y_max - y_min)
        ax.set_ylim(y_min - pad, y_max + pad)

    return fig

def main() -> int:
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
