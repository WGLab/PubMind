import re

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

    variant = re.sub(r"\s+", "", variant)
    variant = variant.replace("→", ">").replace("c.", "").replace("C.", "").replace("g.", "").replace("G.", "")

    patterns = [
        (r"\b(\d+)([ATGC]+)[-/]*[>/]([ATGC]+)\b", lambda m: (m[0], m[1], m[2], None)),
        (r"\b([ATGC]+)(\d+)([ATGC]+)\b", lambda m: (m[1], m[0], m[2], None)),
        (r"\b(\d+)del([ATGC]*)\b", lambda m: (m[0], m[1] if m[1] else "-", "del", None)),
        (r"\bdel(\d+)([ATGC]*)\b", lambda m: (m[0], m[1] if m[1] else "-", "del", None)),
        (r"\b(\d+)_\d+del([ATGC]*)\b", lambda m: (m[0], m[1] if m[1] else "-", "del", None)),
        (r"\b(\d+)ins([ATGC]+)\b", lambda m: (m[0], "ins", m[1], None)),
        (r"\bins(\d+)([ATGC]+)\b", lambda m: (m[0], "ins", m[1], None)),
        (r"\b(\d+)dup([ATGC]+)\b", lambda m: (m[0], "dup", m[1], None)),
        (r"\bdup(\d+)([ATGC]+)\b", lambda m: (m[0], "dup", m[1], None)),
        (r"\b(rs\d+)([ATGC]+)[-/]*[>/]([ATGC]+)\b", lambda m: (None, m[1], m[2], m[0])),
        (r"\b(rs\d+)\b", lambda m: (None, None, None, m[0]))
    ]

    for pattern, handler in patterns:
        match = re.search(pattern, variant, re.IGNORECASE)
        if match:
            return handler(match.groups())

    return None, None, None, None


# Evaluate
###############Protein############
lst_of_test = ['p.Lys1129Ter','p.Lys1129*','P.p.Lys1129Ter',
               'p.(lys121arg)','Lys121arg','Lys121Arg',
               'A21*','T32Ter','T342sp','K32stop',
               'A12K','K32N',
               'A21del','Lys32del','K3DEL',
               '21insA','32insLys','3INSK',
               'A21fs','Lys32fs',
               'A21fs*','Lys32fs*','A21fs*21','Lys32fs*31',
               '(Pro588Argfs*18)'
              ]

lst_of_test_results = ['Lys1129Ter','Lys1129Ter','Lys1129Ter',
               'lys121arg','Lys121arg','Lys121Arg',
               'Ala21Ter','Thr32Ter','Thr342Ter','Lys32Ter',
               'Ala12Lys','Lys32Asn',
               'Ala21del','Lys32del','Lys3del',
               'ins21Ala','ins32Lys','ins3Lys',
               'Ala21fs','Lys32fs',
               'Ala21fs*','Lys32fs*','Ala21fs*21','Lys32fs*31',
               'Pro588Argfs*18'
              ]

results = []
for i, test_str in enumerate(lst_of_test):
    parsed = parse_protein_variant(test_str)
    try:
        if not parsed[0] or not parsed[1] or not parsed[2]:
            result_str = "None"
        elif parsed[1] == "ins":
            result_str = f"ins{parsed[0]}{parsed[2]}"
        else:
            result_str = f"{parsed[1]}{parsed[0]}{parsed[2]}"
        results.append((test_str, result_str, lst_of_test_results[i],result_str == lst_of_test_results[i]))
    except Exception as e:
        results.append((test_str, f"ERROR: {e}", False))

if len([re for re in results if not re[3]])==0:
    print("Pass all tests!")
else:
    print("Failed cases:")
    print([re for re in results if not re[3]])
    
###########cDNA###########
# Example test cases and expected results
test_variants = [
    "c.123A>G", "c.(123A>G)","G123A", "123delA", "del123T", 
    "123_456del", "123insA", "ins123C",
    "123dupT", "dup123G", "rs123", "RS1029524",
    "rs456C>T", "g.456T>G", "  321  A  >  T"
]
expected_outputs = [
    ("123", "A", "G", None), ("123", "A", "G", None), ("123", "G", "A", None), ("123", "A", "del", None), ("123", "T", "del", None),
    ("123", "-", "del", None), ("123", "ins", "A", None), ("123", "ins", "C", None),
    ("123", "dup", "T", None), ("123", "dup", "G", None), (None, None, None, "rs123"),(None, None, None, "RS1029524"),
    (None, "C", "T", "rs456"), ("456", "T", "G", None), ("321", "A", "T", None)
]

# Evaluate function
results = []
for var, expected in zip(test_variants, expected_outputs):
    try:
        parsed = parse_dna_variant(var)
        results.append((var, parsed, expected,parsed == expected))
    except Exception as e:
        results.append((var, f"ERROR: {e}", False))

if len([re for re in results if not re[3]])==0:
    print("Pass all tests!")
else:
    print("Failed cases:")
    print([re for re in results if not re[3]])  
    

    
    
#########Extracting LLM output##########
#load
df = pd.read_csv('../PubVarDB/normlization/1-initial_quality_filter/run_20250425.genename.tsv',sep='\t',
                 header=0,dtype=str,low_memory=False)

### patho filter
df["pathogenicity"] = df_all["pathogenicity"].replace("-", "unknown")
# refine patho classification
valid_patho = ['pathogenic', 'likely pathogenic', 'conflicting', 'likely benign', 'benign', 'unknown']
df = df[df['pathogenicity'].isin(valid_patho)]


# Apply functions to extract normalized columns
df[['dna_pos', 'dna_from', 'dna_to', 'RSID']] = df['DNA mutation'].apply(lambda x: pd.Series(parse_dna_variant(x)))
df[['aa_pos', 'aa_from', 'aa_to']] = df['protein mutation'].apply(lambda x: pd.Series(parse_protein_variant(x)))
df['dna_change']=df['dna_from']+df['dna_pos']+df['dna_to']
df['aa_change']=df['aa_from']+df['aa_pos']+df['aa_to']


#save
df.to_csv('../PubVarDB/normlization/2-variant_regex_extract/run_20250425.genename.patho_filter.variant_norm.tsv',sep='\t',index=False)