# PubVarDB
PubVarDB is a human genetic variant database generated using large language models (LLMs) applied to millions of PubMed abstracts and PubMed Central full-text articles.
![PubVarDB Overview copy](https://github.com/user-attachments/assets/0dac7052-6995-4b9c-952e-e3876627ca50)


## Remaining Tasks (ranked in Priority)
- [X] Evalaute and finalize BERT-filtering (get larger train dataset with ~10000 abstracts, try BioMedBERT and PubMedBERT, compare final LLM output)
- [ ] Do a seperate LLM run using the paragraph that filtered out by old distillBERT but kept by BioMedBERT
- [ ] Finalize the DB (incorperate protein change -> genome coordinates), run ANNOVAR and check model performance
- [ ] Build Web Access (using SQL + Phen2Gene template)
- [ ] Write a Paper
- [ ] Try using different prompt to generate seperated DBs for CNV, gene fushion, long indel, SV.
