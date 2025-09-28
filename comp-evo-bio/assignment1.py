#!/usr/bin/env python3
import argparse
import re
import sys

def parse_args():
    parser = argparse.ArgumentParser(description="Scale branch lengths in a Newick tree.")
    parser.add_argument("-i", "--input", required=True, help="Input Newick file (e.g., input.tre)")
    parser.add_argument("-s", "--scale", required=True, type=float, help="Scale factor (float)")
    parser.add_argument("-o", "--output", required=True, help="Output Newick file (e.g., scaled.tre)")
    return parser.parse_args()

def scale_newick(newick_str, scale):
    # Regex to find numbers after a colon (branch lengths)
    def scale_match(match):
        length = float(match.group(1))
        return f":{length * scale:.8f}"  # keep precision with 8 decimals
    return re.sub(r":([0-9.]+)", scale_match, newick_str)

def main():
    args = parse_args()

    # Step 1: Read input file
    try:
        with open(args.input, "r") as f:
            newick_str = f.read().strip()
    except FileNotFoundError:
        print(f"Error: Cannot read file {args.input}", file=sys.stderr)
        sys.exit(1)

    # Step 2: Scale branch lengths
    try:
        scaled_str = scale_newick(newick_str, args.scale)
    except Exception as e:
        print(f"Error while scaling: {e}", file=sys.stderr)
        sys.exit(1)

    # Step 3: Write output
    try:
        with open(args.output, "w") as f:
            f.write(scaled_str + "\n")
    except Exception as e:
        print(f"Error: Cannot write to {args.output}: {e}", file=sys.stderr)
        sys.exit(1)

    # Step 4: Print original + scaled trees
    print("Original tree:")
    print(newick_str)
    print("\nScaled tree:")
    print(scaled_str)

if __name__ == "__main__":
    main()