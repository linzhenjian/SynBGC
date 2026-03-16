#!/usr/bin/env python3
import sys

if len(sys.argv) != 3:
    print(f"Usage: {sys.argv[0]} <prefix_file> <full_file>", file=sys.stderr)
    sys.exit(1)

prefix_file, full_file = sys.argv[1], sys.argv[2]

# Load prefixes into a set
prefixes = set()
with open(prefix_file) as f:
    for line in f:
        p = line.strip()
        if p:
            prefixes.add(p)

# Match full names that start with any prefix
with open(full_file) as f:
    for line in f:
        name = line.strip()
        for prefix in prefixes:
            if name.startswith(prefix):
                print(name)
                break  # matched one prefix, no need to check others

