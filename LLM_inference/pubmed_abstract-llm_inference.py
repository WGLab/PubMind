from bs4 import BeautifulSoup
import pandas as pd
from transformers import AutoTokenizer
import glob
from vllm import LLM, SamplingParams
import re
import os
import sys

filename_list_path=sys.argv[1] #for multiple slurm jobs submission

def generate_prompt(tmvar_corpus, few_shot_examples):
    
    prompts= []
    for paragraph in tmvar_corpus:
        prompt=f'{few_shot_examples}\nParagraph: \"{paragraph}\"'
        messages = [
            {'role': 'system', 'content': 
             '''You are an expert in genomics and bioinformatics. Extract the variant-related information from scientific publication paragraphs.'''},
            {'role': 'user', 'content': prompt},
        ]

        # Apply chat template
        templated_prompt = tokenizer.apply_chat_template(
            messages,
            add_generation_prompt=True,  # Adds necessary tokens for generation
            tokenize=False               # Returns a string
        )
        
        templated_prompt=templated_prompt.replace('Cutting Knowledge Date: December 2023\nToday Date: 26 Jul 2024\n\n','')
        prompts.append(templated_prompt)
    
    return prompts

# Function to convert model output to a DataFrame
def inference_to_df(output_test, pmc_ids):
    rows = []
    for i, x in enumerate(output_test):
        pmc_id = pmc_ids[i]  # Assign the correct PMC ID
        #sec = sec_for_inputs[i]
        for line in x.outputs[0].text.strip().splitlines():
            if "##" in line: # allow ## not at the beginning
                parts = line.split("##", 1)[-1].strip().split("||") 
                rows.append([pmc_id] + parts)
    
    df_cot = pd.DataFrame(rows)
    df_cot = df_cot.iloc[:,:8]
    #if df_cot.shape[1] == 7:  # Ensure the columns match the expected number
    df_cot.columns = ['PMID', 'gene', 'DNA mutation', 'protein mutation', 'disease', 'phenotype', 'LLM reasoning', 'pathogenicity']
    df_cot = df_cot.drop_duplicates()
    return df_cot


### 20250429
snv_fewshot_prompts = """
Extract the information of gene, DNA mutation (or RSID), protein mutation, related disease, related phenotype, reason for pathogenicity, and pathogenicity from this research paper paragraph. Answer exactly in the format of ##gene||DNA mutation||protein mutation||disease||phenotype||pathogenicity reasoning||pathogenicity. Make sure to choose the pathogenicity only within these catogories [pathogenic, likely pathogenic, conflicting, likely benign, benign, unknown].

Example 1:
Paragraph: "Paragraph: We used zebra fish model to test some candidate variants for orofacial cleft (OFC). For ESRP2, variants R250Q and R667C rescued the molecular splicing of Arhgef11 in the Py2T assay, (Fig. 3E, Table 1)."
Output: 
##ESRP2||-||R250Q||orofacial cleft (OFC)||-||rescued the molecular splicing of Arhgef11 in the Py2T assay||benign
##ESRP2||-||R667C||orofacial cleft (OFC)||-||rescued the molecular splicing of Arhgef11 in the Py2T assay||benign

Example 2:
Paragraph: "The mutation c.123A>T in BRCA1 leads to a truncated protein p.Lys41*, associated with breast cancer. Functional assays showed loss of DNA repair function, indicating pathogenicity."
Output: 
##BRCA1||c.123A>T||p.Lys41*||breast cancer||-||loss of DNA repair function||pathogenic

Example 3:
Paragraph: "A variant in CFTR, p.Phe508del, causes cystic fibrosis. Studies confirmed defective chloride channels, confirming pathogenicity."
Output: 
##CFTR||-||p.Phe508del||cystic fibrosis||defective chloride channels||defective chloride channels||pathogenic
"""

###20250714
gene_fusion_fewshot_prompt="""Extract gene fusion information from the paragraph in the format:  
##gene fusion||driver gene||partner gene||domain affected||related disease||experiment or functional result||pathogenicity

- The **gene fusion** should be in the format **GENE1::GENE2** as mentioned in the text.
- Identify the **driver gene** — this is typically the gene contributing the active domain (e.g., kinase, transcription factor) or the 3' gene in oncogenic fusions. Use context clues (e.g., "drives", "activates", "encodes") to determine.
- The **partner gene** is the fusion partner that contributes less functional consequence or is upstream/regulatory.

Example 1:  
"The EML4-ALK fusion occurs in lung cancer and activates the ALK kinase domain, leading to cell proliferation in vitro."  
##EML4::ALK||EML4||ALK||kinase domain||lung cancer||cell proliferation in vitro||pathogenic

Example 2:  
"A rare fusion between TMPRSS2 and ERG was found in prostate cancer, with no functional assay performed."  
##TMPRSS2::ERG||TMPRSS2||ERG||-||prostate cancer||not tested||unknown

Example 3:  
"We identified a novel RUNX1-RUNX1T1 fusion in acute myeloid leukemia, disrupting the DNA-binding domain of RUNX1."  
##RUNX1::RUNX1T1||RUNX1||RUNX1T1||DNA-binding domain||acute myeloid leukemia||disrupts DNA binding||pathogenic

Now extract from the paragraph below:

"""

###20250714
cnv_fewshot_prompt="""
You are given a biomedical paragraph. Extract any copy number variants (CNVs) mentioned, including deletions or duplications. For each CNV, extract the following structured fields:

##CNV_type||chromosome_region||gene(s)||genomic_coordinates||disease||experiment or functional result||pathogenicity

- CNV_type: "deletion" or "duplication"
- chromosome_region: cytogenetic band (e.g., 22q11.2) if available
- gene(s): affected gene(s), separated by commas if multiple
- genomic_coordinates: if available (e.g., chr22:18000000-21000000)
- disease: associated disease name, or "unknown"
- experiment or functional result: brief summary of experimental or clinical findings
- pathogenicity: must be one of ["pathogenic", "likely pathogenic", "benign", "likely benign", "unknown", "conflicting"]

Example 1:
"The 22q11.2 deletion syndrome is caused by a ~3Mb deletion on chromosome 22, encompassing the TBX1 gene, and is associated with DiGeorge syndrome."
Answer:
##deletion||22q11.2||TBX1||chr22:18000000-21000000||DiGeorge syndrome||3Mb deletion spanning TBX1||pathogenic

Example 2:
"A microduplication of 16p11.2 involving the SH2B1 gene was observed in individuals with early-onset obesity."
Answer:
##duplication||16p11.2||SH2B1||unknown||early-onset obesity||Reported in multiple obesity patients||likely pathogenic

Example 3:
"An inherited 1q21.1 deletion encompassing the GJA5 gene was identified in a family with variable phenotypes, but with no consistent clinical presentation."
Answer:
##deletion||1q21.1||GJA5||unknown||unknown||Detected in unaffected individuals and patients||conflicting

Example 4:
"A 2q13 deletion involving NPHP1 was found in both patients with kidney disease and unaffected family members."
Answer:
##deletion||2q13||NPHP1||unknown||kidney disease||Lack of consistent phenotype among carriers||likely benign


Now extract the CNV(s) from the following paragraph:
"""

###20250714
sv_fewshot_prompt = """
You are given a biomedical paragraph. Extract any **structural variants (SVs)** mentioned, including deletions, duplications, inversions, insertions, and translocations. For each SV, extract the following structured fields:

##SV_type||gene(s)||chromosome_region||genomic_coordinates||related disease||experimental or clinical evidence||pathogenicity

- SV_type: one of ["deletion", "duplication", "inversion", "insertion", "translocation"]
- gene(s): affected gene(s), separated by commas if multiple
- chromosome_region: cytogenetic region (e.g., 22q11.2), or "unknown"
- genomic_coordinates: format "chr:start-end" if available, or "unknown"
- related disease: disease associated with the SV, or "unknown"
- experimental or clinical evidence: brief summary of supporting data
- pathogenicity: one of ["pathogenic", "likely pathogenic", "benign", "likely benign", "unknown", "conflicting"]

Example 1:
"A pathogenic inversion disrupting the F8 gene on Xq28 is a common cause of severe hemophilia A."
Answer:
##inversion||F8||Xq28||chrX:154180000-154300000||hemophilia A||Breakpoints identified in patients with severe phenotype||pathogenic

Example 2:
"A translocation between chromosomes 9 and 22 creates the BCR-ABL1 fusion gene, known as the Philadelphia chromosome, found in chronic myeloid leukemia."
Answer:
##translocation||BCR,ABL1||9q34,22q11||chr9:133729000-chr22:23632600||chronic myeloid leukemia||Philadelphia chromosome detected by FISH||pathogenic

Example 3:
"A 5q13.3 duplication involving the SMN1 and SMN2 genes has been reported in some individuals without symptoms."
Answer:
##duplication||SMN1,SMN2||5q13.3||chr5:70000000-70200000||unknown||Identified in asymptomatic carriers during screening||benign

Example 4:
"An insertion of ~1.5 kb in intron 1 of the FGFR2 gene has been observed in individuals with Apert syndrome."
Answer:
##insertion||FGFR2||10q26||chr10:123800000||Apert syndrome||Insertion disrupts regulatory elements in FGFR2 intron 1||likely pathogenic

Now extract the SV(s) from the following paragraph:
"""

#########!!!!!!!!!!#########
#########define the prompt you want to run here########
chosen_prompt=snv_fewshot_prompts
#########!!!!!!!!!!#########

# Initialize the model
model = LLM('Llama-3.3-70B-Instruct',
    tensor_parallel_size=2, # Change this to match the number of GPUs
    max_model_len=5000,
    gpu_memory_utilization=0.95,
    swap_space=0
)

tokenizer = AutoTokenizer.from_pretrained("Llama-3.3-70B-Instruct")

# Base directory for csv/tsv batched input
batch_pmc_path = '/mnt/isilon/wang_lab/pengwang/projects/LLM/PubVarDB/BERT_filtered_input_pubmed'

llm_run_batch_lookup_path = '/mnt/isilon/wang_lab/pengwang/projects/LLM/PubVarDB/llm_run_batch_lookup_files'
# Read the file names from 'run1_1-90' or other system provided path
with open(os.path.join(llm_run_batch_lookup_path, filename_list_path), 'r') as file:
    filenames = file.readlines()

for filename in filenames:
    filename = filename.strip()
    if filename.endswith(".xml"): 
        file_path = os.path.join(batch_pmc_path, filename)
        print(f"Loading: {filename}")
        
        if os.path.exists(file_path) and os.path.getsize(file_path) > 2:
            batch=pd.read_csv(file_path,sep='\t')
            
        else:
            print(f"File {file_path} is empty or does not exist.")
            continue # move to the next file
            
        all_inputs=list(batch['Abstract'])
        final_inputs = generate_prompt(all_inputs,chosen_prompt)
        pmc_ids_for_inputs=list(batch['PMID'])
        #sec_for_inputs=list(batch['section'])

        print("Running inference ...")
        # Run inference on all inputs
        sampling_params = SamplingParams(temperature=0.6, top_p=0.3, max_tokens=2000) #1000, increase to 2000 in case
        output = model.generate(final_inputs, sampling_params)

        # Convert inference results to a single DataFrame
        df_cot = inference_to_df(output, pmc_ids_for_inputs)

        # Save the results to a TSV file
        output_file = f"/mnt/isilon/wang_lab/pengwang/projects/LLM/PubVarDB/LLM_output/LLM_output_pubmed/{filename.replace('.xml', '.out.tsv')}"
        df_cot.to_csv(output_file, header=True, index=False, sep='\t')

        print(f"Saved results to {output_file}")
