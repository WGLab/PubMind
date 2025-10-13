# <img width="180" height="70" alt="pubmind_logo_v1" src="https://github.com/user-attachments/assets/077c2f07-33a5-4cd6-9680-ec38e8ea4b21" />


PubMind is a large language model (LLM)-assisted framework for <ins>**Pub**</ins>lication <ins>**M**</ins>utation and <ins>**in**</ins>formation <ins>**D**</ins>iscovery, designed to extract variant–disease–pathogenicity relationships directly from biomedical literature.

<img width="1095" alt="image" src="https://github.com/user-attachments/assets/fa4717dd-07e6-48d4-8a50-64ecd478b807" />

PubMind is an AI-driven framework that uses large language models (LLMs) to extract genetic variant–disease–pathogenicity associations directly from biomedical literature. It combines fine-tuned BERT models for input filtering with instruction-tuned LLMs for extracting variant, disease, and functional evidence, covering SNVs, CNVs, SVs, and gene fusions. Extracted variants are normalized to genomic and transcript coordinates and stored in PubMind-DB, a web-accessible knowledgebase. Applied to >41M PubMed abstracts and >5M PMC full texts, PubMind-DB contains ~0.7M consolidated unique variants with rich annotations, of which only ~10% overlap with ClinVar—yet >80% of those show concordant pathogenicity labels, including full agreement for four-star expert-reviewed variants. PubMind provides a scalable, generalizable, and open-source framework that transforms unstructured text into structured genomic knowledge, supporting variant interpretation and precision medicine.

## Prerequisite

Please refer to `requirements.txt` for required packages.

## Run PubMind

Please refer to `run_PubMind.ipynb` for how to use PubMind. All inputs and outputs during this example PubMind run are in the `example` folder.

PubMind frameworkds includes the following modules:
1. Filtering Module (finetuned BERT model)
   - Wangwpi/PubMind_finetuned_BERT (Hugging Face)
2. Inference Module (instruction-tuned LLM)
   - meta-llama/Llama-3.3-70B-Instruct (Hugging Face)
3. Normalization Module
   - Quality filter (gene name, pathogenicity)
   - Variant parser (cDNA, protein, RSID)
   - Map to transcript
   - Map to genome cooridnates
   - MONDO Disease name
   - HPO term
   
## PubMind-DB

PubMind-DB could be accessed here: https://pubmind.wglab.org/

## License

PubMind is freely available for academic use. For license details, please refer to [this page](LICENSE.md).
