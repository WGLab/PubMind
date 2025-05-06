import pandas as pd
import re
import glob

### human gene name
from intervaltree import IntervalTree
import pyensembl
from Bio.Seq import Seq
from Bio import pairwise2, SeqIO
from pyensembl import EnsemblRelease
data = EnsemblRelease(release=111, species='homo_sapiens')
gene_names = data.gene_names()
gene_names = [gene for gene in gene_names if gene !='']

### pubmed
folder_path = "../PubVarDB/LLM_output/run_20250425/LLM_output_pubmed/"  # Update this with your actual folder path
tsv_files = glob.glob(folder_path + "*.tsv")
df_list = [pd.read_csv(file, sep="\t") for file in tsv_files]
df = pd.concat(df_list, ignore_index=True)
df=df.rename(columns={'disease':'phenotype','phenotype':'disease'})

df = df[df['gene'].isin(gene_names)]
df=df[(df['DNA mutation']!='-') | (df['protein mutation']!='-')]

### pmcid-pmid transversion
pmcid_pmid = pd.read_csv('../pubmed/pmid_pmcid.lookup',sep='\t',header=None,dtypes=str)
pmcid_pmid.columns=['PMCID','PMID']
pmcid_pmid['PMID']=pmcid_pmid['PMID'].apply(lambda x: str(x).replace('PMID:',''))

df_merge=df.merge(pmcid_pmid,on='PMID',how='left')
df_merge=df_merge[['PMCID','PMID', 'gene', 'DNA mutation', 'protein mutation', 'phenotype',
       'disease', 'LLM reasoning', 'pathogenicity']]

### PMC
folder_path = "../PubVarDB/LLM_output/run_20250425/LLM_output_all_pmc/"  # Update this with your actual folder path
tsv_files = glob.glob(folder_path + "*.tsv")

df2_list = [pd.read_csv(file, sep="\t") for file in tsv_files]
df2 = pd.concat(df2_list, ignore_index=True)

df2=df2.rename(columns={'disease':'phenotype','phenotype':'disease'})
df2 = df2[df2['gene'].isin(gene_names)]
df2=df2[(df2['DNA mutation']!='-') | (df2['protein mutation']!='-')]

df2_merge = df2.merge(pmcid_pmid,on='PMCID',how='left')[['PMCID', 'PMID','section', 'gene', 'DNA mutation', 'protein mutation',
       'phenotype','disease', 'LLM reasoning', 'pathogenicity']]

### combine PMC and pubmed
df_all = pd.concat([df2_merge,df_merge]).drop_duplicates()
df_all['pathogenicity'] = df_all['pathogenicity'].apply(lambda x: str(x).strip())

df_all.to_csv('../PubVarDB/normlization/1-initial_quality_filter/run_20250425.genename.tsv',sep='\t',index=False)