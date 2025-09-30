# PubMind
PubMind is a large language model (LLM)-assisted framework for <ins>**Pub**</ins>lication <ins>**M**</ins>utation and <ins>**in**</ins>formation <ins>**D**</ins>iscovery, designed to extract variant–disease–pathogenicity relationships directly from biomedical literature.

<img width="1095" alt="image" src="https://github.com/user-attachments/assets/fa4717dd-07e6-48d4-8a50-64ecd478b807" />

PubMind is for academic use. For commercial use please contact CHOP office of technology transfer.

# Running PubMind

## Prerequisite
Please refer to `requirements.txt` for necessary python packages.

## Use PubMind
Please refer to `run_PubMind.ipynb` for the use of PubMind. All example input and output are in the `example` folder.

PubMind includes the following steps:
1. Filtering Module (finetuned BERT model)
2. Inference Module (instruction-tuned LLM)
3. Normalization Module
   - Quality filter (gene name, pathogenicity)
   - Variant parser (cDNA, protein, RSID)
   - Map to transcript
   - Map to genome cooridnates
   - MONDO Disease name
   - HPO term
   
