#!/usr/bin/env python3
"""
Plot the ONT raw signal for a single read from a .pod5 file.

Usage:
    python plot_pod5_signal.py <file.pod5> <read_id_uuid>

Dependencies:
    pip install pod5 numpy matplotlib
"""

import argparse
import sys
import uuid
import numpy as np
import matplotlib.pyplot as plt

try:
    import pod5 as p5
except Exception as e:
    print("Error: the 'pod5' package is required. Install with 'pip install pod5'.", file=sys.stderr)
    raise

def main() -> int:
    parser = argparse.ArgumentParser(
        description="Plot signal vs time for a specific read_id from a POD5 file."
    )
    parser.add_argument("pod5_path", help="Path to the .pod5 file")
    parser.add_argument("read_id", help="UUID read_id to plot (e.g. 0000173c-bf67-...)")
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

            # Pretty plot
            plt.figure(figsize=(10, 3))
            plt.plot(t, signal_pa, linewidth=0.5)
            plt.xlabel("Time (s)")
            plt.ylabel("Current (pA)")

            # Optional metadata in the title if available
            ch = getattr(getattr(read, "pore", None), "channel", None)
            title_bits = [f"read {read_id}", f"{int(sample_rate)} Hz", f"{n} samples"]
            if ch is not None:
                title_bits.insert(1, f"channel {ch}")
            plt.title(" | ".join(title_bits))

            plt.tight_layout()
            plt.show()

    except FileNotFoundError:
        print(f"Error: file not found: {args.pod5_path}", file=sys.stderr)
        return 1
    except Exception as e:
        print(f"Unhandled error: {e}", file=sys.stderr)
        return 99

    return 0

if __name__ == "__main__":
    raise SystemExit(main())