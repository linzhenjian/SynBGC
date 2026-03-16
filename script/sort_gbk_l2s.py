#!/usr/bin/env python3

import sys
from Bio import SeqIO

inp, outp = sys.argv[1], sys.argv[2]

records = list(SeqIO.parse(inp, "genbank"))
records.sort(key=lambda r: len(r.seq), reverse=True)

SeqIO.write(records, outp, "genbank")

