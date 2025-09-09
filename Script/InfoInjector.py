
import os
import argparse
import pod5
import uuid
from dataclasses import replace

# -----------------------------------------------------------------------------
# 1. Define Constant to be written into run info per sample
# -----------------------------------------------------------------------------
FLOWCELL_CODE = "FLO-MIN114"
SAMPLE_ID = "SIMULATED"

# -----------------------------------------------------------------------------
# 2. Argument parsing with argparse
# -----------------------------------------------------------------------------
def parse_args():
    """
    Parse command-line arguments for input POD5 and output POD5 paths.
    """
    parser = argparse.ArgumentParser(
        description=(
            "Read a simulated POD5 file, fix missing run_info and channel_info fields "
            "for Dorado compatibility, and write out a new POD5."
        )
    )
    parser.add_argument("inputfile", type=str)
    parser.add_argument("outputfile", type=str)
    return parser.parse_args()

# -----------------------------------------------------------------------------
# 3. Main logic: open Reader on the input POD5, Writer on the output POD5
# -----------------------------------------------------------------------------
def main():
    args = parse_args()
    input_pod5  = args.inputfile
    if not os.path.isfile(input_pod5):
        raise argparse.ArgumentError(f"The file '{input_pod5}' does not exist.") # type: ignore
    output_pod5 = args.outputfile

    # Open the input POD5 for reading
    with pod5.Reader(input_pod5) as reader, pod5.Writer(output_pod5) as writer:
        for record in reader:
            read = record.to_read()  # Convert immutable record → mutable Read

            # -----------------------------------------------------------------
            # Populate missing run_info fields so Dorado can basecall properly
            # -----------------------------------------------------------------
            tracking = dict(getattr(read.run_info, "tracking_id", {}) or {})
            tracking["run_id"] = str(uuid.uuid4())
            tracking.setdefault("exp_start_time", "0.0")
            ctx = dict(getattr(read.run_info, "context_tags", {}) or {})
            ctx["sample_frequency"] = "5000"
            updated_info = replace(
                read.run_info,
                flow_cell_product_code=FLOWCELL_CODE,   # :contentReference[oaicite:8]{index=8}
                sample_id=SAMPLE_ID,
                tracking_id=tracking,
                context_tags=ctx
            )
            read.run_info = updated_info

            # Write the edited read to the output POD5
            writer.add_read(read)
    print(f"✓ Finished processing. Output POD5 written to: {output_pod5}")

if __name__ == "__main__":
    main()