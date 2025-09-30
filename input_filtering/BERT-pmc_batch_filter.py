# prepare a dataset for BERT using LLM
from bs4 import BeautifulSoup
import pandas as pd
from transformers import AutoConfig
import glob
import re
import os
import time
import torch
import sys
from transformers import AutoModelForSequenceClassification, AutoTokenizer

#process each run based on sys input (should need 27 runs in total for 263 files)
n_run=int(sys.argv[1])

# Run prediction using BERT
def predict(text):
    inputs = tokenizer(text, truncation=True, padding="max_length", max_length=512, return_tensors="pt")
    
    # Move input tensors to the same device as the model
    inputs = {key: val.to(device) for key, val in inputs.items()}
    
    # Disable gradient calculation for inference
    with torch.no_grad():
        outputs = model(**inputs)
    
    logits = outputs.logits
    prediction = torch.argmax(logits, dim=-1).item()
    return prediction

# Run BERT model to decide keep a paragraph or not
def keep_paragraph_bert(para):
    pred = predict(para)
    if pred==1:
        return True
    else:
        return False

# Function to extract all paragraphs along with their sections and abstract
def get_full_text_paragraphs(i_pmc):
    directory_path = f"/mnt/isilon/wang_lab/shared/pubmed_central_additional/pmc_nxml/{i_pmc}/"
    nxml_files = glob.glob(f"{directory_path}*.nxml")

    if nxml_files:
        with open(nxml_files[0], "r", encoding="utf-8", errors="replace") as file:
            content = file.read()
        soup = BeautifulSoup(content, "xml")
        
        
        # Extract paragraphs
        paragraphs = []
        
        abstract_section = soup.find("abstract")
        if abstract_section:
            abstract_text = abstract_section.get_text(separator=" ", strip=True)
            if keep_paragraph_bert(abstract_text):
                paragraphs.append({"section": "abstract", "paragraph": abstract_text})
            
        for sec in soup.find_all("sec"):
            section_title = sec.find("title").get_text(strip=True) if sec.find("title") else "No Title"

            # Skip irrelevant sections
            section_not_interest = ["method", "acknowledgment", "reference"]
            if any(sec_name in section_title.lower() for sec_name in section_not_interest):
                continue
            
            valid_sections = ["abstract", "introduction", "method", "result", "discussion", "conclusion"
                             "methods", "results", "discussions"]
            if not any(sec == section_title.lower() for sec in valid_sections):
                    section_title = "Others"
            for p in sec.find_all("p"):
                paragraph_text = p.get_text(separator=" ", strip=True)
                if keep_paragraph_bert(paragraph_text):
                    paragraphs.append({"section": section_title, "paragraph": paragraph_text})

        # Convert to DataFrame
        df = pd.DataFrame(paragraphs,columns=['section','paragraph'])
        
        unknown_section_mask = df["section"] == "Others"

        # Drop duplicates only in "Others"
        df_sort = pd.concat([
            df[~unknown_section_mask],  
            df[unknown_section_mask].drop_duplicates(subset=["paragraph"]) # make sure "Others" rank at the bottom
        ])

        df_sort = df_sort.drop_duplicates(subset=["paragraph"], keep="first").reset_index(drop=True)

        return df_sort

    return pd.DataFrame(columns=["section", "paragraph"])

# Process a batch of PMC IDs from a file
def process_batch_file(batch_file, output_directory):
    # Read PMC IDs from the batch file
    with open(batch_file, "r") as f:
        pmc_ids = f.read().splitlines()

    df_all=pd.DataFrame(columns=['PMC_ID','section','paragraph'])
    
    count=0
    for pmc_id in pmc_ids:
        count += 1
        df = get_full_text_paragraphs(pmc_id)
        df['PMC_ID']=pmc_id
        
        df_all = pd.concat([df_all,df])
        
        if count%1000 == 0: #total of 10k papers, print 10 times for each 1k finish
            print('.',end='')

    # Save the processed batch to a new CSV file in the output directory
    batch_name = os.path.basename(batch_file).replace(".txt", ".tsv")
    output_file = os.path.join(output_directory, batch_name)
    df_all.to_csv(output_file, index=False, sep='\t')
    print(f"\nProcessed batch saved to {output_file}")

print(f'Run Number: {n_run}')

# Input and output directories
batch_directory = "/mnt/isilon/wang_lab/pengwang/projects/LLM/pubmed_central_additional/batch_pmc_ids/"  # Adjust the path to your batch files
output_directory = "/mnt/isilon/wang_lab/pengwang/projects/LLM/pubmed_central_additional/batch_BERT_filtered_input/"  # Path to save filtered results

# Create output directory if it doesn't exist
os.makedirs(output_directory, exist_ok=True)
batch_files = set(glob.glob(os.path.join(batch_directory, "*.txt")))

#creating the batch for each run (parallel running, each use 1 GPU)
prefix_run = [str(i)+'-' for i in list(range((n_run-1)*10+1,n_run*10+1))]


# Extract and filter based on prefix
batch_files_run = [f for f in batch_files if os.path.basename(f).startswith(tuple(prefix_run))]

# Load the BERT model
model = AutoModelForSequenceClassification.from_pretrained("/mnt/isilon/wang_lab/pengwang/projects/LLM/BERT/saved_model")

# Load the tokenizer
tokenizer = AutoTokenizer.from_pretrained("/mnt/isilon/wang_lab/pengwang/projects/LLM/BERT/saved_model")

# Move to GPU if available
device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
model.to(device)

print("Model successfully loaded!")


for batch_file in batch_files_run:
    t0=time.time()
    print(f'Processing: {os.path.basename(batch_file)}')
    process_batch_file(batch_file, output_directory)
    t1=time.time()
    print(f"Finished time: {t1-t0}.")
    print('='*15)