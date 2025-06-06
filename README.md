# PubMind
PubMind is a large language model (LLM)-assisted framework for <ins>**Pub**</ins>lication <ins>**M**</ins>utation and <ins>**in**</ins>formation <ins>**D**</ins>iscovery, designed to extract variant–disease–pathogenicity relationships directly from biomedical literature.

<img width="1095" alt="image" src="https://github.com/user-attachments/assets/fa4717dd-07e6-48d4-8a50-64ecd478b807" />



## Remaining Tasks (ranked in Priority)
- [X] Evalaute and finalize BERT-filtering (get larger train dataset with ~10000 abstracts, try BioMedBERT and PubMedBERT, compare final LLM output)
- [X] Finalize the DB (incorperate protein change -> genome coordinates), run ANNOVAR and check model performance
- [X] Build Web Access (using SQL + Phen2Gene template)
- [X] Write a Paper
- [ ] Try using different prompt to generate seperated DBs for CNV, gene fushion, long indel, SV.
