#!/usr/bin/env python3
"""
filter_gbk_preserve_origin_line.py

• Keeps only contigs that have >5 valid gene–CDS pairs
• Trims feature qualifiers (gene + optional translation)
• Removes “_” from contig names in LOCUS / ACCESSION / VERSION
  while preserving the DEFINITION line and the DNA sequence
"""

import sys
import os
import re
from Bio import SeqIO

# ────────────────────────────────────────────────────────────────────────
def filter_gbk_preserve_origin_line(input_gbk: str, output_gbk: str) -> None:
    cleaned_records = []   # SeqRecord objects that will be written

    # ── Step 1: parse GenBank records one by one ────────────────────────
    for record in SeqIO.parse(input_gbk, "genbank"):
        gene_feature_by_id = {}   # gene_id ➜ gene feature
        cds_kept_for_gene = {}    # gene_id ➜ *one* CDS feature

        # Pass 1 – collect gene/CDS and pair them by shared id
        for feat in record.features:
            if feat.type not in {"gene", "CDS"}:
                continue

            # prefer /gene; fallback to /locus_tag
            gene_id = feat.qualifiers.get(
                "gene", feat.qualifiers.get("locus_tag", [None])
            )[0]
            if not gene_id:
                continue

            if feat.type == "gene":
                gene_feature_by_id[gene_id] = feat
                cds_kept_for_gene.setdefault(gene_id, None)
            else:  # CDS
                if cds_kept_for_gene.get(gene_id) is None:
                    cds_kept_for_gene[gene_id] = feat

        # Pass 2 – keep only well-paired gene/CDS features
        new_features = []
        for gid, gene_feat in gene_feature_by_id.items():
            cds_feat = cds_kept_for_gene.get(gid)
            if cds_feat is None:
                continue

            # derive a gene name
            gene_name = gene_feat.qualifiers.get("gene", [None])[0]
            if not gene_name:
                for key in ("locus_tag", "note", "product"):
                    val = gene_feat.qualifiers.get(key)
                    if val:
                        gene_name = val[0]
                        break
            if not gene_name:
                continue  # cannot name this pair ➜ skip

            # tidy qualifiers
            gene_feat.qualifiers = {"gene": [gene_name]}
            cds_qualifiers = {"gene": [gene_name]}
            if "translation" in cds_feat.qualifiers:
                cds_qualifiers["translation"] = cds_feat.qualifiers["translation"]
            cds_feat.qualifiers = cds_qualifiers

            new_features.extend((gene_feat, cds_feat))

        # skip contigs with ≤5 pairs (each pair = 2 features)
        if len(new_features) // 2 <= 5:
            continue

        record.features = new_features     # keep DNA; do NOT strip record.seq

        # ── Remove underscores from contig name everywhere it matters ───
        new_id = record.id.replace("_", "")
        record.id = new_id
        record.name = new_id      # affects LOCUS
        # DO NOT touch record.description → DEFINITION remains intact

        if "accessions" in record.annotations:
            record.annotations["accessions"] = [
                acc.replace("_", "") for acc in record.annotations["accessions"]
            ]
        if "sequence_version" in record.annotations:
            ver = record.annotations["sequence_version"]
            if isinstance(ver, str):
                ver = ver.replace("_", "")
            record.annotations["sequence_version"] = ver

        cleaned_records.append(record)

    # safety check
    if not cleaned_records:
        sys.exit("No contigs with >5 genes having CDS found; nothing written.")

    # ── Step 2: write temporary GenBank file ────────────────────────────
    tmp_file = output_gbk + ".tmp"
    SeqIO.write(cleaned_records, tmp_file, "genbank")

    # ── Step 3: patch LOCUS lines (remove underscores) & write final file
    with open(tmp_file) as fin, open(output_gbk, "w") as fout:
        for line in fin:
            if line.startswith("LOCUS"):
                m = re.match(r"^(LOCUS\s+)(\S+)(\s+.*)$", line)
                if m:
                    prefix, locus_id, suffix = m.groups()
                    locus_id = locus_id.replace("_", "")
                    line = f"{prefix}{locus_id}{suffix}\n"
            fout.write(line)

    os.remove(tmp_file)
    print(
        f"[DONE] GenBank written to: {output_gbk} "
        f"({len(cleaned_records)} contigs kept, DNA preserved)"
    )

# ── CLI ─────────────────────────────────────────────────────────────────
if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python filter_gbk_preserve_origin_line.py input.gbk output.gbk")
        sys.exit(1)
    filter_gbk_preserve_origin_line(sys.argv[1], sys.argv[2])

