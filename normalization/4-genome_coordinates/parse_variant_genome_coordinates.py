import pandas as pd
import ast
from Bio.Seq import Seq
from Bio.Data import CodonTable

# 3-letter to 1-letter amino acid dictionary
aa_dict = {
    'A': 'Ala', 'R': 'Arg', 'N': 'Asn', 'D': 'Asp', 'C': 'Cys',
    'Q': 'Gln', 'E': 'Glu', 'G': 'Gly', 'H': 'His', 'I': 'Ile',
    'L': 'Leu', 'K': 'Lys', 'M': 'Met', 'F': 'Phe', 'P': 'Pro',
    'S': 'Ser', 'T': 'Thr', 'W': 'Trp', 'Y': 'Tyr', 'V': 'Val',
    '*': 'Ter'
}
aa_3to1 = {v: k for k, v in aa_dict.items()}

# Amino acid to RNA codon table
standard_table = CodonTable.unambiguous_rna_by_name["Standard"]
aa_to_codons = {}
for codon, aa in standard_table.forward_table.items():
    aa_to_codons.setdefault(aa, []).append(codon)
for codon in standard_table.stop_codons:
    aa_to_codons.setdefault("*", []).append(codon)


def get_cdna_variant(tx, dna_pos, dna_from, dna_to):
    if 'ref_cDNA' not in tx or 'genomic_coord' not in tx:
        return None

    # Compare reference base at cDNA level
    if str(tx['ref_cDNA']).upper() != str(dna_from).upper():
        return None

    # genomic_coord is a single position (int), not a list
    pos = tx['genomic_coord']
    ref = str(dna_from).upper()
    alt = str(dna_to).upper()

    # Adjust for strand
    if tx['strand'] == '-':
        ref = str(Seq(ref).complement())
        alt = str(Seq(alt).complement())

    return {
        'chr': tx['chromosome'],
        'pos': pos,
        'ref': ref,
        'alt': alt,
        'transcript_id': tx['transcript_id']
    }

    
    
def get_protein_variant(tx, aa_pos, aa_from, aa_to):
    if 'ref_aa' not in tx or 'genomic_coords' not in tx or 'gnomic_codon' not in tx:
        return []

    # Convert 3-letter to 1-letter amino acid
    aa_from_1 = aa_3to1.get(str(aa_from).capitalize())
    aa_to_1 = aa_3to1.get(str(aa_to).capitalize())
    if not aa_from_1 or not aa_to_1:
        return []

    if tx['ref_aa'].upper() != aa_from_1:
        return []

    ref_codon = tx['gnomic_codon'].upper()
    strand = tx['strand']
    coords = tx['genomic_coords']
    variants = []

    for alt_rna_codon in aa_to_codons.get(aa_to_1, []):
        # Convert RNA to DNA
        alt_dna_codon = alt_rna_codon.replace("U", "T")

        if strand == '+':
            alt_genomic_codon = alt_dna_codon
            start_coord = coords[0]
        else:
            alt_genomic_codon = str(Seq(alt_dna_codon).reverse_complement())
            start_coord = coords[-1]

        if ref_codon != alt_genomic_codon:
            variants.append({
                'chr': tx['chromosome'],
                'pos': start_coord,  # always the leftmost genomic position
                'ref': ref_codon,
                'alt': alt_genomic_codon,
                'transcript_id': tx['transcript_id'],
                'aa_from': aa_from,
                'aa_to': aa_to
            })

    return variants


def parse_variant_from_transcript(row, mane_tx_dict):
    gene = row['gene']
    transcripts = ast.literal_eval(row['genomic_coord_result'])

    # Extract variant information
    dna_pos = row.get('dna_pos')
    dna_from = row.get('dna_from')
    dna_to   = row.get('dna_to')
    aa_pos   = row.get('aa_pos')
    aa_from  = row.get('aa_from')
    aa_to    = row.get('aa_to')

    mane_tx = mane_tx_dict.get(gene)
    results = []

    for tx in transcripts:
        is_mane = (tx['transcript_id'] == mane_tx)

        if 'ref_aa' in tx and pd.notnull(aa_pos) and pd.notnull(aa_from) and pd.notnull(aa_to):
            variants = get_protein_variant(tx, aa_pos, aa_from, aa_to)
        elif 'ref_cDNA' in tx and pd.notnull(dna_pos) and pd.notnull(dna_from) and pd.notnull(dna_to):
            variant = get_cdna_variant(tx, dna_pos, dna_from, dna_to)
            variants = [variant] if variant else []
        else:
            continue

        for var in variants:
            if var:
                var_clean = {
                    'chr': var['chr'],
                    'pos': var['pos'],
                    'ref': var['ref'],
                    'alt': var['alt'],
                    'MANE_transcript_used': str(is_mane)
                }
                results.append(var_clean)

    # Deduplicate by chr, pos, ref, alt
    seen = set()
    unique_results = []
    for r in results:
        key = (r['chr'], r['pos'], r['ref'], r['alt'])
        if key not in seen:
            seen.add(key)
            unique_results.append(r)

    return unique_results if unique_results else None

