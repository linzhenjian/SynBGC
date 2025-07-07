#!/usr/bin/env python3

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature, FeatureLocation, CompoundLocation
import argparse
import re

def find_bounds(features):
    starts, ends = [], []
    for feat in features:
        if feat.type in {"gene", "CDS"}:
            loc = feat.location
            if isinstance(loc, CompoundLocation):
                for part in loc.parts:
                    starts.append(int(part.start))
                    ends.append(int(part.end))
            else:
                starts.append(int(loc.start))
                ends.append(int(loc.end))
    if not starts or not ends:
        return None
    return min(starts), max(ends)

def filter_qualifiers(feat):
    new_qual = {}
    if feat.type == "gene":
        if "gene" in feat.qualifiers:
            new_qual["gene"] = feat.qualifiers["gene"]
    elif feat.type == "CDS":
        for key in ("gene", "translation", "note"):
            if key in feat.qualifiers:
                new_qual[key] = feat.qualifiers[key]
    return new_qual

def shift_location(location, offset):
    if isinstance(location, CompoundLocation):
        new_parts = [
            FeatureLocation(int(part.start) - offset, int(part.end) - offset, strand=part.strand)
            for part in location.parts
        ]
        return CompoundLocation(new_parts, operator="join")
    else:
        return FeatureLocation(
            int(location.start) - offset,
            int(location.end) - offset,
            strand=location.strand
        )

def trim_record(record, prefix=""):
    bounds = find_bounds(record.features)
    if not bounds:
        return None
    start, end = bounds

    new_seq = record.seq[start:end]
    new_features = []

    for feat in record.features:
        if feat.type not in {"gene", "CDS"}:
            continue

        loc = feat.location
        if isinstance(loc, CompoundLocation):
            if any(part.start < start or part.end > end for part in loc.parts):
                continue
        else:
            if loc.start < start or loc.end > end:
                continue

        new_loc = shift_location(loc, offset=start)
        new_feat = SeqFeature(location=new_loc, type=feat.type, qualifiers=filter_qualifiers(feat))
        new_features.append(new_feat)

    # Extract {species}_contig format and apply prefix
    match = re.match(r"(\{.*?\}_)?(.+)", record.id)
    species_part = match.group(1) or ""
    base_name = match.group(2)
    new_id = f"{species_part}{prefix}{base_name}"

    new_record = SeqRecord(
        seq=new_seq,
        id=new_id,
        name=new_id,
        description=record.description,
        annotations=record.annotations.copy(),
        features=new_features
    )

    # Update annotations to reflect new LOCUS name
    new_record.annotations["accessions"] = [new_id]
    new_record.annotations["source"] = record.annotations.get("source", "")
    new_record.annotations["organism"] = record.annotations.get("organism", "")

    return new_record

def trim_gbk(input_gbk, output_gbk, prefix):
    trimmed = []
    for record in SeqIO.parse(input_gbk, "genbank"):
        print(f"Processing {record.id}")
        cleaned = trim_record(record, prefix=prefix)
        if cleaned:
            trimmed.append(cleaned)
        else:
            print(f"  Skipped {record.id}, no CDS/gene in bounds")

    if trimmed:
        SeqIO.write(trimmed, output_gbk, "genbank")
        print(f"✅ Trimmed GenBank written to: {output_gbk}")
    else:
        print("❌ No valid records to write.")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Trim GenBank file: keep gene/CDS with selected qualifiers and rename LOCUS.")
    parser.add_argument("input_gbk", help="Input GenBank file")
    parser.add_argument("output_gbk", help="Output trimmed GenBank file")
    parser.add_argument("--prefix", default="", help="Prefix to add to contig name in LOCUS line")
    args = parser.parse_args()
    trim_gbk(args.input_gbk, args.output_gbk, args.prefix)

