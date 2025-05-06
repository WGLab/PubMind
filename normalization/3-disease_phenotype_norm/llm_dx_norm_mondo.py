import json
import os
import pandas as pd
from sentence_transformers import SentenceTransformer, util
from difflib import SequenceMatcher
import numpy as np
import sys

mondo_df = pd.read_csv('../../../MONDO/MONDO_id_name.tsv',sep='\t',header=0)

job_id=int(sys.argv[1])

mondo_dx_lst = list(mondo_df['MONDO name'])

emb_model = SentenceTransformer('neuml/pubmedbert-base-embeddings') ## Change your embedding model here
features_embeddings = emb_model.encode(mondo_dx_lst, convert_to_tensor=True)

def text_similarity(a, b):
    return SequenceMatcher(None, a, b).ratio()  # Ratio between 0 and 1

def find_similarity_emb(entity, graph_features_nodes,features_embeddings, emb_model = emb_model, threshold = 0.9):
    # Compute embeddings for the entity and graph_nodes
    entity = entity.lower()
    #if graph_features_nodes:
    #    features_embeddings = emb_model.encode(graph_features_nodes, convert_to_tensor=True)

    # features_embeddings = features_embeddings.to(self.device)
    entity_embedding = emb_model.encode(entity, convert_to_tensor=True)
    # Compute cosine similarities
    similarities = util.cos_sim(entity_embedding, features_embeddings).squeeze(0)
    # Filter graph_nodes by similarity threshold
    similar_terms =  [(node,sim) for node, sim in zip(graph_features_nodes, similarities) if sim >= threshold]#,
    if len(similar_terms) == 0:
        similar_terms =  [(node,sim) for node, sim in zip(graph_features_nodes, similarities) if sim >= threshold-0.05]#,

    #similar_terms = sorted(similar_terms, key=lambda x: x[1], reverse=True)
    # Convert all tensors to float in one step (avoids repeated `.item()` calls)
    semantic_scores = np.array([sim.cpu().item() for _, sim in similar_terms])  # Move to CPU and convert to NumPy array
    text_similarities = np.array([text_similarity(entity, term) for term, _ in similar_terms])

    # Weighted combination of scores
    alpha = 0.75  # Adjust to control balance
    final_scores = alpha * semantic_scores + (1 - alpha) * text_similarities  # Vectorized computation

    # Sort indices based on final scores in descending order
    sorted_indices = np.argsort(final_scores)[::-1]

    # Extract sorted results
    sorted_terms = [(similar_terms[i][0], semantic_scores[i], text_similarities[i], final_scores[i]) for i in sorted_indices]
    if len(sorted_terms) > 0:
        return sorted_terms[0] #only return the highest term
    else:
        return ()
    

# load llm dx output
llm_out = pd.read_csv('../2-variant_regex_extract/run_20250425.genename.patho_filter.variant_norm.tsv',
                     sep='\t',dtype='str')['disease'].unique()

llm_dx_lst = list(llm_out)

# Split data into 4 chunks
num_splits = 4
chunk_size = len(llm_dx_lst) // num_splits
start_idx = job_id * chunk_size
end_idx = (job_id + 1) * chunk_size if job_id < num_splits - 1 else len(llm_dx_lst)
llm_dx_lst_job = llm_dx_lst[start_idx:end_idx]

print('Run ID:', job_id)
print('Total LLM disease to be normalized:', len(llm_dx_lst_job))


# perform normalization
dx_llm_mondo_dict = {}
counter=1
for dx in llm_dx_lst_job:
    
    if counter % 5==0:
        print(counter,end=' ')
        
    counter+=1
    
    dx_llm_mondo_dict[dx]=find_similarity_emb(dx,mondo_dx_lst,features_embeddings,emb_model,0)
    
print()
print('Finish normalization, saving the result...')

df_lst=[[k,v[0],v[1],v[2],v[3]] if len(v)!=0 else [k,None,np.nan,np.nan,np.nan] for k,v in dx_llm_mondo_dict.items() ]
df_norm = pd.DataFrame(df_lst,columns=['disease','MONDO name','disease pbmbert sim','disease string sim','disease combine sim'])
df_norm_merge = df_norm.merge(mondo_df,on='MONDO name',how='left')
df_norm_merge.to_csv(f'llm_mondo_norm_{job_id}.tsv', sep='\t',index=False,header=True)
print(f'Disease normalization saved to llm_mondo_norm_{job_id}.tsv')