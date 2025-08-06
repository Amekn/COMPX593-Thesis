#!/usr/bin/env python3
"""
compute_mean_qscore_template_percentage.py

Script to compute the average of the 'mean_qscore_template' column as a percentage for all sequences
in a Nanopore summary TSV file.

Usage:
    python compute_mean_qscore_template_percentage.py <report_file> [--convert]

Arguments:
    report_file: Path to the TSV report file with a header that includes 'mean_qscore_template'.
Options:
    -c, --convert   Convert Phred Q-scores into accuracy percentages via:
                    accuracy = (1 - 10**(-Q/10)) * 100

If --convert is not set, the script simply reports the mean Q‑score.
"""
import argparse
import sys
import pandas as pd

def parse_args():
    p = argparse.ArgumentParser(
        description="Compute the average mean_qscore_template as a percentage from a TSV report."
    )
    p.add_argument(
        'report_file',
        help='TSV file path containing mean_qscore_template'
    )
    p.add_argument(
        '-c', '--convert',
        action='store_true',
        help='Convert Phred Q-scores to accuracy percentages'
    )
    return p.parse_args()


def main():
    args = parse_args()
    try:
        df = pd.read_csv(args.report_file, sep='\t')
    except Exception as e:
        print(f"Error reading '{args.report_file}': {e}", file=sys.stderr)
        sys.exit(1)

    col_name = 'mean_qscore_template'
    if col_name not in df.columns:
        print(f"Missing required column: '{col_name}'", file=sys.stderr)
        sys.exit(1)

    # Convert to numeric and drop non-numeric rows
    qvals = pd.to_numeric(df[col_name], errors='coerce').dropna()
    if qvals.empty:
        print(f"No valid numeric entries in '{col_name}'", file=sys.stderr)
        sys.exit(1)

    if args.convert:
        # Compute accuracy from the **mean** Q‑score
        mean_q = qvals.mean()
        mean_pct = (1 - 10 ** (-mean_q / 10)) * 100
        print(f"Mean accuracy percentage based on mean Q‑score ({mean_q:.2f}): {mean_pct:.2f}%")
    else:
        # Report the mean Q‑score directly
        mean_q = qvals.mean()
        print(f"Mean Q‑score of '{col_name}': {mean_q:.2f}")

if __name__ == '__main__':
    main()
