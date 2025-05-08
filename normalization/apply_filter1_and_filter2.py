import pandas as pd
import re
import glob
import sys
import os
from pyensembl import EnsemblRelease

input_path = sys.argv[1]  # First argument: input folder path
output_path=input_path.replace('.out.tsv','.out.qc.tsv')

#2-variant_regex_extract
# Dictionary to convert one-letter to three-letter amino acid codes
aa_dict = {
    'A': 'Ala', 'R': 'Arg', 'N': 'Asn', 'D': 'Asp', 'C': 'Cys',
    'Q': 'Gln', 'E': 'Glu', 'G': 'Gly', 'H': 'His', 'I': 'Ile',
    'L': 'Leu', 'K': 'Lys', 'M': 'Met', 'F': 'Phe', 'P': 'Pro',
    'S': 'Ser', 'T': 'Thr', 'W': 'Trp', 'Y': 'Tyr', 'V': 'Val',
    '*': 'Ter'  # Stop codon
    
}

def parse_protein_variant(variant):
    """Parse protein variant and return (position, from_aa, to_aa)"""
    if not isinstance(variant, str) or variant in ["-", ""]:
        return None, None, None

    # Normalize input
    variant = re.sub(r"\s+", "", variant)
    variant = variant.replace("p.", "").replace("P.", "")
    variant = variant.replace("Ter", "*").replace("sp", "*").replace("stop", "*")
    variant = variant.replace("DEL", "del").replace("Del", "del").replace("INS", "ins").replace("Ins", "ins")

    patterns = [
        # 3-letter substitution + frameshift + *
        (r"([A-Za-z]{3})(\d+)([A-Za-z]{3})fs\*(\d+)", lambda m: (m[1], aa_dict.get(m[0], m[0]), f"{m[2]}fs*{m[3]}")),
        # 3-letter substitution
        (r"\b([A-Za-z]{3})(\d+)([A-Za-z]{3})\b", lambda m: (m[1], m[0], m[2])),
        # 1-letter substitution
        (r"\b([A-Z])(\d+)([A-Z])\b", lambda m: (m[1], aa_dict.get(m[0], m[0]), aa_dict.get(m[2], m[2]))),
        # frameshift + stop
        (r"([A-Za-z]{3})(\d+)fs\*(\d+)", lambda m: (m[1], aa_dict.get(m[0], m[0]), f"fs*{m[2]}")),
        (r"([A-Z])(\d+)fs\*(\d+)", lambda m: (m[1], aa_dict.get(m[0], m[0]), f"fs*{m[2]}")),
        # frameshift + no stop
        (r"([A-Za-z]{3})(\d+)fs\*", lambda m: (m[1], aa_dict.get(m[0], m[0]), "fs*")),
        (r"([A-Z])(\d+)fs\*", lambda m: (m[1], aa_dict.get(m[0], m[0]), "fs*")),
        # simple frameshift
        (r"\b([A-Za-z]{3})(\d+)fs\b", lambda m: (m[1], aa_dict.get(m[0], m[0]), "fs")),
        (r"\b([A-Z])(\d+)fs\b", lambda m: (m[1], aa_dict.get(m[0], m[0]), "fs")),
        # deletion
        (r"\b([A-Za-z]{3})(\d+)del\b", lambda m: (m[1], aa_dict.get(m[0], m[0]), "del")),
        (r"\b([A-Z])(\d+)del\b", lambda m: (m[1], aa_dict.get(m[0], m[0]), "del")),
        # termination
        (r"([A-Za-z]{3})(\d+)\*", lambda m: (m[1], aa_dict.get(m[0], m[0]), "Ter")),
        (r"([A-Z])(\d+)\*", lambda m: (m[1], aa_dict.get(m[0], m[0]), "Ter")),
        # insertion
        (r"\b(\d+)ins([A-Z])\b", lambda m: (m[0], "ins", aa_dict.get(m[1], m[1]))),
        (r"\b(\d+)ins([A-Za-z]{3})\b", lambda m: (m[0], "ins", aa_dict.get(m[1], m[1])))
    ]

    for pattern, handler in patterns:
        match = re.search(pattern, variant)
        if match:
            return handler(match.groups())

    return None, None, None

# Concise version of parse_dna_variant
def parse_dna_variant(variant):
    """Return (pos, from, to, rsid) from a DNA/cDNA mutation string."""
    if not isinstance(variant, str) or variant.strip() in ["-", ""]:
        return None, None, None, None

    variant = re.sub(r"\s+", "", variant) #removes all whitespace characters
    variant = variant.replace("→", ">").replace("c.", "").replace("C.", "").replace("g.", "").replace("G.", "")

    patterns = [
        (r"\b([\d\+]+)([ATGC]+)[-/]*[>/]([ATGC]+)\b", lambda m: (m[0], m[1], m[2], None)),
        (r"\b(\d+)([ATGC]+)[-/]*[>/]([ATGC]+)\b", lambda m: (m[0], m[1], m[2], None)),
        (r"\b([ATGC]+)(\d+)([ATGC]+)\b", lambda m: (m[1], m[0], m[2], None)),
        (r"\b(\d+)del([ATGC]*)\b", lambda m: (m[0], m[1] if m[1] else "-", "del", None)),
        (r"\bdel(\d+)([ATGC]*)\b", lambda m: (m[0], m[1] if m[1] else "-", "del", None)),
        (r"\b(\d+_\d+)del([ATGC]*)\b", lambda m: (m[0], m[1] if m[1] else "-", "del", None)),
        (r"\b([\d_]+)ins([ATGC]+)\b", lambda m: (m[0], "ins", m[1], None)),
        (r"\bins([\d_]+)([ATGC]+)\b", lambda m: (m[0], "ins", m[1], None)),
        (r"\b([\d_]+)dup([ATGC]+)\b", lambda m: (m[0], "dup", m[1], None)),
        (r"\bdup([\d_]+)([ATGC]+)\b", lambda m: (m[0], "dup", m[1], None)),
        (r"\b(rs\d+)([ATGC]+)[-/]*[>/]([ATGC]+)\b", lambda m: (None, m[1], m[2], m[0])),
        (r"\b(rs\d+)\b", lambda m: (None, None, None, m[0]))
    ]

    for pattern, handler in patterns:
        match = re.search(pattern, variant, re.IGNORECASE)
        if match:
            return handler(match.groups())

    return None, None, None, None

data = EnsemblRelease(release=111, species='homo_sapiens')
gene_names = data.gene_names()
gene_names = [gene for gene in gene_names if gene !='']

#1-initial_quality_filter
df = pd.read_csv(input_path, sep="\t")
df = df[df['gene'].isin(gene_names)]
df=df[(df['DNA mutation']!='-') | (df['protein mutation']!='-')]
df.to_csv(output_path,sep='\t')
print(f'Finished filter-1: save output at {output_path}')


#2-variant_regex_extract
### patho filter
df["pathogenicity"] = df["pathogenicity"].replace("-", "unknown")
# refine patho classification
valid_patho = ['pathogenic', 'likely pathogenic', 'conflicting', 'likely benign', 'benign', 'unknown']
df = df[df['pathogenicity'].isin(valid_patho)]

# Apply functions to extract normalized columns
df[['dna_pos', 'dna_from', 'dna_to', 'RSID']] = df['DNA mutation'].apply(lambda x: pd.Series(parse_dna_variant(x)))
df[['aa_pos', 'aa_from', 'aa_to']] = df['protein mutation'].apply(lambda x: pd.Series(parse_protein_variant(x)))
df['dna_change']=df['dna_from']+df['dna_pos']+df['dna_to']
df['aa_change']=df['aa_from']+df['aa_pos']+df['aa_to']

#save
output_path2 = output_path.replace('.out.qc.tsv','.out.qc.var_regex.tsv')
df.to_csv(output_path2,sep='\t',index=False)
print(f'Finished filter-2: saved to {output_path2}')