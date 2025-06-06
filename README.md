# PubMind
PubMind is an LLM-based framework for **Pub**lication **M**utation and **in**formation **D**iscovery directly using PubMed abstracts and central full texts.

<img width="468" alt="image" src="https://github.com/user-attachments/assets/174bf110-efcf-4ceb-93f3-25f7129c56ea" />



## Remaining Tasks (ranked in Priority)
- [X] Evalaute and finalize BERT-filtering (get larger train dataset with ~10000 abstracts, try BioMedBERT and PubMedBERT, compare final LLM output)
- [X] Finalize the DB (incorperate protein change -> genome coordinates), run ANNOVAR and check model performance
- [X] Build Web Access (using SQL + Phen2Gene template)
- [X] Write a Paper
- [ ] Try using different prompt to generate seperated DBs for CNV, gene fushion, long indel, SV.
