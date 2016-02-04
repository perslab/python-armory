#!/usr/bin/python

import numpy as np
import pandas as pd
import pdb

# Function to map from Mm symbol to Hs Ensembl gene identifiers
def get_mapping(df_mm2mm,df_mm2hs,mm_symbol):
	pdb.set_trace()
	#mm_symbol = dge_gene_id.split(":")[-1]
	if mm_symbol in df_mm2mm.index:

		# Discard many-to-many mappings (e.g. Pou6f1, which maps to ENSMUSG00000009739 and ENSMUSG00000098598)
		matches_mm = df_mm2mm.ix[mm_symbol,'Ensembl Gene ID']
		if isinstance(matches_mm, str) and matches_mm in df_mm2hs.index: # isinstance check fails if several matches for mm_symbol

			# Keep one-to-many mappings, by selecting most idential human gene
			matches_hs = df_mm2hs.ix[matches_mm,:]
			if len(matches_hs.shape) > 1:
				row_index = np.argmax(matches_hs.ix[:,"% Identity with respect to Human gene"].tolist()) 
				return matches_hs.ix[row_index,'Human Ensembl Gene ID'] 
			else:
				return matches_hs['Human Ensembl Gene ID'] 
	return "not_found"
	
# Function to convert mouse gene identifiers to Ensembl human identifiers
def to_ensembl(df_mm2mm,df_mm2hs,df):
	rows2drop = []
	mapping = []
	for ix in df.index.tolist():
		hs_ensembl_id = get_mapping(df_mm2mm,df_mm2hs,ix)
		if hs_ensembl_id == "not_found":
			rows2drop.append(ix)
		else:
			mapping.append(hs_ensembl_id)
	df.drop(rows2drop,axis=0,inplace=True)
	df['Ensembl Gene ID'] = pd.Series(mapping, index=df.index)
	return df.groupby('Ensembl Gene ID',sort=False).mean(),rows2drop

