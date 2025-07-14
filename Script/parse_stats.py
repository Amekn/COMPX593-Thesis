#!/usr/bin/env python3
"""
parse_stats.py: Parse a .stats input (from a file or stdin) and display its contents in a formatted table.

Usage:
    # Read from a file:
    python parse_stats.py sup_exp.stats

    # Read from stdin:
    cat sup_exp.stats | python parse_stats.py

Dependencies:
    pip install tabulate
"""
import sys
import argparse
from tabulate import tabulate

def parse_stats(file_obj):
    """
    Read a .stats file-like object and extract metric, value, and optional comment.
    Lines starting with 'SN' are parsed; others are ignored.
    """
    rows = []
    for line in file_obj:
        line = line.strip()
        # Only process lines beginning with 'SN\t'
        if not line.startswith("SN\t"):
            continue
        parts = line.split("\t")
        # parts: ['SN', metric, value, comment?]
        metric = parts[1]
        value = parts[2]
        comment = ""
        if len(parts) > 3:
            # Remove leading '#' if present
            comment = parts[3].lstrip('#').strip()
        rows.append((metric, value, comment))
    return rows

def main():
    parser = argparse.ArgumentParser(
        description="Parse a .stats input (file or stdin) into a formatted table"
    )
    parser.add_argument(
        "stats_file",
        nargs="?",
        type=argparse.FileType('r'),
        default=sys.stdin,
        help="Path to the .stats file; if omitted, reads from stdin"
    )
    args = parser.parse_args()

    # Parse rows from the provided file or stdin
    rows = parse_stats(args.stats_file)
    if not rows:
        print("No valid 'SN' entries found in the input.")
        return

    # Display the table with custom alignment
    print('\nAlignment Statistics:\n')
    print(tabulate(
        rows,
        headers=["Metric", "Value", "Comment"],
        tablefmt="fancy_grid",
        colalign=("left", "right", "left")
    ))

if __name__ == "__main__":
    main()
