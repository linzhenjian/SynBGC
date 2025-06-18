#!/usr/bin/env python3
from Bio import SeqIO
from Bio.Seq import Seq
import sys, os

def filter_gbk_preserve_origin_line(input_gbk, output_gbk):
    # ── Step 1: grab every LOCUS line in the original file ────────────────
    original_locus_lines = [
        line.rstrip() for line in open(input_gbk) if line.startswith("LOCUS")
    ]

    cleaned_records = []     # records to keep
    kept_locus_lines = []    # corresponding LOCUS lines to preserve

    # ── Step 2: parse GenBank records one by one ─────────────────────────
    for rec_idx, record in enumerate(SeqIO.parse(input_gbk, "genbank")):
        gene_feature_by_id = {}     # stores gene features
        cds_kept_for_gene = {}      # stores one CDS per gene

        # ── Pass 1: Identify and pair gene/CDS features ───────────────────
        for feat in record.features:
            if feat.type not in {"gene", "CDS"}:
                continue  # ignore other features

            # Get gene identifier: prefer /gene, fallback to /locus_tag
            gene_id = feat.qualifiers.get("gene", feat.qualifiers.get("locus_tag", [None]))[0]
            if not gene_id:
                continue  # skip if no usable ID

            if feat.type == "gene":
                gene_feature_by_id[gene_id] = feat
                cds_kept_for_gene.setdefault(gene_id, None)

            elif feat.type == "CDS":
                if cds_kept_for_gene.get(gene_id) is None:
                    cds_kept_for_gene[gene_id] = feat  # keep only the first CDS per gene

        # ── Step 3: Validate, clean and keep matching gene-CDS pairs ──────
        new_features = []
        for gid, gene_feat in gene_feature_by_id.items():
            cds_feat = cds_kept_for_gene.get(gid)
            if cds_feat is None:
                continue  # skip gene if no CDS found

            # Try to get /gene qualifiers
            gene_name = gene_feat.qualifiers.get("gene", [None])[0]
            cds_name = cds_feat.qualifiers.get("gene", [None])[0]

            # Try to infer name if missing in either
            if not gene_name or not cds_name:
                shared_keys = set(gene_feat.qualifiers) & set(cds_feat.qualifiers)
                for key in ["locus_tag", "note", "product"]:
                    if key in shared_keys:
                        inferred = gene_feat.qualifiers[key][0]
                        gene_name = cds_name = inferred
                        break

            # Still no gene name? skip this pair
            if not gene_name:
                continue

            # Set /gene tag consistently
            gene_feat.qualifiers["gene"] = [gene_name]
            cds_feat.qualifiers["gene"] = [gene_name]

            # Clean up qualifiers: keep only /gene (and /translation for CDS)
            gene_feat.qualifiers = {"gene": [gene_name]}
            cds_qualifiers = {"gene": [gene_name]}
            if "translation" in cds_feat.qualifiers:
                cds_qualifiers["translation"] = cds_feat.qualifiers["translation"]
            cds_feat.qualifiers = cds_qualifiers

            # Add features to the new list
            new_features.extend((gene_feat, cds_feat))

        # ── Step 4: skip contigs with ≤5 valid gene-CDS pairs ─────────────
        if len(new_features) // 2 <= 5:  # each pair adds 2 features
            continue

        record.features = new_features
        record.seq = Seq("")  # strip DNA sequence
        cleaned_records.append(record)
        kept_locus_lines.append(original_locus_lines[rec_idx])

    # ── Step 5: safety check ─────────────────────────────────────────────
    if not cleaned_records:
        sys.exit("No contigs with >5 genes having CDS found; nothing written.")

    # ── Step 6: write intermediate file ──────────────────────────────────
    tmp_file = output_gbk + ".tmp"
    SeqIO.write(cleaned_records, tmp_file, "genbank")

    # ── Step 7: post-process to restore LOCUS and strip ORIGIN lines ─────
    with open(tmp_file) as fin, open(output_gbk, "w") as fout:
        locus_idx = 0
        in_origin = False
        for line in fin:
            if line.startswith("LOCUS"):
                fout.write(kept_locus_lines[locus_idx] + "\n")
                locus_idx += 1
            elif line.startswith("ORIGIN"):
                fout.write("ORIGIN\n")
                in_origin = True
            elif line.strip() == "//":
                fout.write("//\n")
                in_origin = False
            elif in_origin:
                continue  # skip sequence lines
            else:
                fout.write(line)

    os.remove(tmp_file)
    print(f"[DONE] GenBank written to: {output_gbk} "
          f"({len(cleaned_records)} contigs kept, DNA removed)")

# ── CLI: argument handling ───────────────────────────────────────────────
if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python remove_dna_but_keep_origin.py input.gbk output.gbk")
        sys.exit(1)
    filter_gbk_preserve_origin_line(sys.argv[1], sys.argv[2])

