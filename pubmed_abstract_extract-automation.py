from bs4 import BeautifulSoup
import pandas as pd
from transformers import AutoTokenizer
import glob
from vllm import LLM, SamplingParams
import re
import os
import sys

filename_list_path=sys.argv[1]

def generate_prompt(tmvar_corpus):
    
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

###better
### 20250429
few_shot_examples = """
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

# Initialize the model
model = LLM('/mnt/isilon/wang_lab/pengwang/projects/LLM/Llama-3.3-70B-Instruct',
    tensor_parallel_size=2, # Change this to match the number of GPUs
    max_model_len=5000,
    gpu_memory_utilization=0.95,
    swap_space=0
)

tokenizer = AutoTokenizer.from_pretrained("/mnt/isilon/wang_lab/pengwang/projects/LLM/Llama-3.3-70B-Instruct")

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
        final_inputs = generate_prompt(all_inputs)
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
