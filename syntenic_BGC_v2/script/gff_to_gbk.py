#!/usr/bin/env python3
import argparse
import re
from datetime import datetime
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature, FeatureLocation, CompoundLocation

# Exact codon table from Perl script
CODON_TABLE = {
    'UUU': 'F', 'UUC': 'F', 'UUA': 'L', 'UUG': 'L',
    'CUU': 'L', 'CUC': 'L', 'CUA': 'L', 'CUG': 'L',
    'AUU': 'I', 'AUC': 'I', 'AUA': 'I', 'AUG': 'M',
    'GUU': 'V', 'GUC': 'V', 'GUA': 'V', 'GUG': 'V',
    'UCU': 'S', 'UCC': 'S', 'UCA': 'S', 'UCG': 'S',
    'CCU': 'P', 'CCC': 'P', 'CCA': 'P', 'CCG': 'P',
    'ACU': 'T', 'ACC': 'T', 'ACA': 'T', 'ACG': 'T',
    'GCU': 'A', 'GCC': 'A', 'GCA': 'A', 'GCG': 'A',
    'UAU': 'Y', 'UAC': 'Y', 'UAA': '*', 'UAG': '*',
    'CAU': 'H', 'CAC': 'H', 'CAA': 'Q', 'CAG': 'Q',
    'AAU': 'N', 'AAC': 'N', 'AAA': 'K', 'AAG': 'K',
    'GAU': 'D', 'GAC': 'D', 'GAA': 'E', 'GAG': 'E',
    'UGU': 'C', 'UGC': 'C', 'UGA': '*', 'UGG': 'W',
    'CGU': 'R', 'CGC': 'R', 'CGA': 'R', 'CGG': 'R',
    'AGU': 'S', 'AGC': 'S', 'AGA': 'R', 'AGG': 'R',
    'GGU': 'G', 'GGC': 'G', 'GGA': 'G', 'GGG': 'G'
}

def parse_gff(gff_file):
    genes = {}
    with open(gff_file) as f:
        for line in f:
            if line.startswith('#'):
                continue
            parts = line.strip().split('\t')
            if len(parts) < 9:
                continue
            seqid, source, ftype, start, end, score, strand, phase, attributes = parts

            if ftype not in ['gene', 'CDS']:
                continue

            id_str = re.sub(r'.*(transcript_id|gene_id)\s+"([^"]+).*', r'\2', attributes)
            id_str = re.sub(r'\.t\d+.*', '', id_str)

            if not id_str:
                continue

            if ftype == 'gene':
                genes[id_str] = {
                    'contig': seqid,
                    'strand': strand,
                    'cds': [],
                    'start': int(start),
                    'end': int(end)
                }
            elif ftype == 'CDS':
                if id_str not in genes:
                    genes[id_str] = {
                        'contig': seqid,
                        'strand': strand,
                        'cds': [],
                        'start': int(start),
                        'end': int(end)
                    }
                genes[id_str]['cds'].append({
                    'start': int(start),
                    'end': int(end),
                    'strand': strand
                })
    return genes

def reverse_complement(dna):
    complement = str.maketrans('ATCGatcg', 'TAGCtagc')
    return dna.translate(complement)[::-1]

def translate_cds(cds_features, contig_seq, strand):
    ordered_exons = sorted(cds_features, key=lambda x: x['start'])
    full_dna = ''.join(contig_seq[exon['start'] - 1: exon['end']] for exon in ordered_exons)

    if strand == '-':
        full_dna = reverse_complement(full_dna)

    rna = full_dna.upper().replace('T', 'U')

    protein = ''
    for i in range(0, len(rna), 3):
        codon = rna[i:i+3]
        if len(codon) == 3:
            protein += CODON_TABLE.get(codon, 'X')
    return protein

def create_genbank_entry(contig_id, contig_seq, genes, organism, prefix):
    record = SeqRecord(
        Seq(contig_seq),
        id=contig_id,
        description=f"{organism} annotation",
        annotations={
            'molecule_type': 'DNA',
            'topology': 'linear',
            'date': datetime.now().strftime('%d-%b-%Y').upper(),
            'organism': organism,
            'comment': 'Converted from GFF+FASTA with exact Perl replication'
        }
    )

    for gene_id, gene_data in genes.items():
        cds_features = gene_data.get('cds', [])
        if not cds_features:
            continue

        exon_locations = [
            FeatureLocation(
                exon['start'] - 1,
                exon['end'],
                strand=1 if exon['strand'] == '+' else -1
            )
            for exon in sorted(cds_features, key=lambda x: x['start'])
        ]
        protein_seq = translate_cds(cds_features, contig_seq, gene_data['strand'])

        gene_feature = SeqFeature(
            location=FeatureLocation(
                gene_data['start'] - 1,
                gene_data['end'],
                strand=1 if gene_data['strand'] == '+' else -1
            ),
            type='gene',
            qualifiers={
                'gene': [gene_id],
                'locus_tag': [f"{prefix}_{gene_id}"]
            }
        )

        cds_feature = SeqFeature(
            location=CompoundLocation(exon_locations) if len(exon_locations) > 1 else exon_locations[0],
            type='CDS',
            qualifiers={
                'gene': [gene_id],
                'locus_tag': [f"{prefix}{gene_id}"],
                'product': ['hypothetical protein'],
                'protein_id': [f"{prefix}{gene_id}"],
                'translation': [protein_seq],
                'codon_start': [1]
            }
        )

        record.features.append(gene_feature)
        record.features.append(cds_feature)

    return record

def main():
    parser = argparse.ArgumentParser(description='Exact replication of Perl gff-to-gb converter.')
    parser.add_argument('-g', '--gff', required=True, help='Input GFF3 file')
    parser.add_argument('-f', '--fasta', required=True, help='Input FASTA file')
    parser.add_argument('-o', '--output', required=True, help='Output GenBank file')
    parser.add_argument('--organism', default='ORG', help='Organism name')
    parser.add_argument('--prefix', default='', help='Protein ID prefix')
    args = parser.parse_args()

    contigs = SeqIO.to_dict(SeqIO.parse(args.fasta, 'fasta'))
    genes = parse_gff(args.gff)

    contig_genes = {}
    for gene_id, gene_data in genes.items():
        contig = gene_data['contig']
        if contig not in contig_genes:
            contig_genes[contig] = {}
        contig_genes[contig][gene_id] = gene_data

    records = []
    for contig_id, seq_record in contigs.items():
        if contig_id in contig_genes:
            record = create_genbank_entry(
                contig_id,
                str(seq_record.seq),
                contig_genes[contig_id],
                args.organism,
                args.prefix
            )
            records.append(record)

    with open(args.output, 'w') as out_handle:
        SeqIO.write(records, out_handle, 'genbank')

if __name__ == '__main__':
    main()

