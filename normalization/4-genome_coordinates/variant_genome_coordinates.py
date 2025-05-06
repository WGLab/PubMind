import pandas as pd
from pyensembl import EnsemblRelease
from Bio.Seq import Seq
import numpy as np
import logging
logging.getLogger("pyensembl.sequence_data").setLevel(logging.ERROR)


# Initialize Ensembl reference (GRCh38, release 104)
ensembl = EnsemblRelease(111, species='homo_sapiens')

aa_dict = {
    'A': 'Ala', 'R': 'Arg', 'N': 'Asn', 'D': 'Asp', 'C': 'Cys',
    'Q': 'Gln', 'E': 'Glu', 'G': 'Gly', 'H': 'His', 'I': 'Ile',
    'L': 'Leu', 'K': 'Lys', 'M': 'Met', 'F': 'Phe', 'P': 'Pro',
    'S': 'Ser', 'T': 'Thr', 'W': 'Trp', 'Y': 'Tyr', 'V': 'Val',
    '*': 'Ter'  # Stop codon
}

def get_cds_genomic_positions(tx):
    ranges = tx.coding_sequence_position_ranges
    genomic_positions = []
    for start, end in ranges:
        if tx.strand == "+":
            genomic_positions.extend(range(start, end + 1))
        else:
            # Reverse direction for minus strand
            genomic_positions.extend(range(end, start - 1, -1))
    return genomic_positions


def protein_to_genomic(gene_name, aa_pos, aa_ref=None):
    try:
        #print(len(ensembl.genes_by_name(gene_name)))
        genes = ensembl.genes_by_name(gene_name)
    except ValueError:
        #print(f'gene {gene_name} not found')
        return None  # Gene not found
    
    all_match_tx = []
    for gene in genes:
        for tx in gene.transcripts:
            if tx.is_protein_coding:
                try:
                    cds = tx.coding_sequence
                    aa_seq = Seq(cds).translate()

                    if aa_pos <= len(aa_seq):
                        aa_from_aa_seq = aa_dict.get(aa_seq[aa_pos - 1])
                        if aa_ref:
                            
                            #check if the user provided aa same as ref aa based on pos
                            if aa_ref!=aa_from_aa_seq:
                                #print('ref aa not match')
                                continue
                                
                        codon_start = (aa_pos - 1) * 3

                        # Get full CDS genomic position map
                        cds_genome = get_cds_genomic_positions(tx)

                        # get genome coodirnate based on codon_start
                        codon_coords = cds_genome[codon_start:codon_start+3]
                        # get cDNA codon
                        codon = cds[codon_start:codon_start + 3]
                        # get genome ref
                        if tx.strand=='+':
                            genome_codon = codon
                        else:
                            genome_codon = str(Seq(codon).reverse_complement())

                        all_match_tx.append({
                            "transcript_id": tx.transcript_id,
                            "chromosome": tx.contig,
                            "strand": tx.strand,
                            "genomic_coords": codon_coords,
                            "ref_aa": aa_seq[aa_pos - 1],
                            "codon_cDNA": codon,
                            'gnomic_codon': genome_codon
                        })
                        
                    #else:
                        #print('ref tx is too short')
                except Exception:
                    continue
    if len(all_match_tx) == 0:
        return None
    return all_match_tx


def cdna_to_genomic(gene_name, cdna_pos, ref_nt=None):
    try:
        genes = ensembl.genes_by_name(gene_name)
    except ValueError:
        #print(f'Gene {gene_name} not found')
        return None

    matches = []

    for gene in genes:
        for tx in gene.transcripts:
            if tx.is_protein_coding:
                try:
                    # Get spliced transcript (cDNA sequence)
                    cdna_seq = tx.coding_sequence
                    if cdna_pos > len(cdna_seq):
                        continue
                    
                    nt_from_tx = cdna_seq[cdna_pos - 1]  # 1-based to 0-based
                    
                    # Optional: check ref nucleotide match
                    if ref_nt:
                        if nt_from_tx.upper() != ref_nt.upper():
                            continue

                    # Use spliced offset → genomic mapping
                    cds_genome = get_cds_genomic_positions(tx)
                    
                    genome_coord = cds_genome[cdna_pos - 1]
                    
                    # get genome ref
                    if tx.strand=='+':
                        genome_ref = nt_from_tx
                    else:
                        genome_ref = str(Seq(nt_from_tx).complement())
                        
                    matches.append({
                        "transcript_id": tx.transcript_id,
                        "chromosome": tx.contig,
                        "strand": tx.strand,
                        "genomic_coord": genome_coord,
                        "ref_cDNA": nt_from_tx,
                        "genomic_ref": genome_ref
                    })
                except Exception as e:
                    continue
    if len(matches) == 0:
        return None
    return matches

def annotate_row(row):
    gene = row['gene']

    if gene=='-' or gene=='nan' or gene==np.nan:
        return None
    
    # Use cDNA info first
    if pd.notnull(row.get('dna_pos')) and pd.notnull(row.get('dna_from')):
        try:
            cdna_coord = cdna_to_genomic(gene, int(row['dna_pos']), row['dna_from'])
            if cdna_coord!=None:
                return cdna_coord
           
        except Exception as e:
            return f'cdna_error: {e}'
        
    # Use protein info if available
    elif pd.notnull(row.get('aa_pos')) and pd.notnull(row.get('aa_from')):
        try:
            return protein_to_genomic(gene, int(row['aa_pos']), row['aa_from'])
        except Exception as e:
            return f'protein_error: {e}'
    
    
    
    return None



cdna_to_genomic('PDGFRB',1681,'C')
cdna_to_genomic('PDGFRB',1681)
cdna_to_genomic('ZSWIM7',201+1, 'G')
protein_to_genomic('PDGFRB',561,'Arg')
protein_to_genomic('CTNND1',439,'Arg')
protein_to_genomic('CTNND1',439)