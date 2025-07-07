# split_fasta.py
import sys
from Bio import SeqIO
from pathlib import Path

fasta = sys.argv[1]
num_chunks = int(sys.argv[2])

records = list(SeqIO.parse(fasta, "fasta"))
chunk_size = (len(records) + num_chunks - 1) // num_chunks

for i in range(num_chunks):
    chunk = records[i * chunk_size:(i + 1) * chunk_size]
    if not chunk:
        continue
    out_path = Path(f"{Path(fasta).stem}_part{i+1}.fa")
    SeqIO.write(chunk, out_path, "fasta")
#python split_fasta.py <input_fasta_file> <number_of_chunks>
