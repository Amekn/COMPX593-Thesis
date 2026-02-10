#!/usr/bin/env python3
"""
Count number of signals per read in a .pod5 file.
Requires: pip install pod5
"""
import sys
from pod5.reader import Reader

def main(pod5_file):
    with Reader(pod5_file) as reader:
        for read in reader.reads():
            read_id = read.read_id  # UUID of the read
            signal_len = len(read.signal)
            print(f"{read_id}\t{signal_len}")

if __name__ == "__main__":
    if len(sys.argv) != 2:
        print(f"Usage: {sys.argv[0]} <file.pod5>")
        sys.exit(1)
    main(sys.argv[1])
