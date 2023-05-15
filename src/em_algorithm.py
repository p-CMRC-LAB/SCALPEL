

import pandas as pd
import numpy as np
import argparse
import vaex as vx
from IPython.display import display
from joblib import Parallel, delayed


#Argument Parser
#**************
parser = argparse.ArgumentParser(description='compute EM algorithm...')
parser.add_argument('cell', type=str, help='path of cell file')
parser.add_argument('output_path', metavar='Efip', type=str, help='path of probability transcripts file')
args = parser.parse_args()


#Functions
#*********

def posterior(probs, tr_estimates):
	'''
	Function to impute transcript relative estimate abundances
	'''
	numerator = (probs * tr_estimates)
	return(numerator.div(numerator.sum(axis=1).values[0]))

def em_algorithm(tab, max_iteration=50):
	'''
	Function to comupte em algorithm on transcript table
	'''
	#number of trs
	tr_size = len(tab.transcript_name.unique())

	#get initial estimated abundances starter
	estimated_abundances = np.round([1/tr_size] * tr_size, 2)
	buffer = estimated_abundances
	#pivot_table
	tab = tab.pivot_table(index='umi',columns='transcript_name',values='frag_prob_weighted', fill_value=0)
	#Iterate until convergence (EM standard)
	for i in range(max_iteration):
		#get posterior probs
		estimated_abundances = np.round((tab.groupby('umi', group_keys=False).apply(lambda x: posterior(x, estimated_abundances))).mean(axis=0),2)
		#keep track if different... else stop
		if np.array_equal(estimated_abundances, buffer):
			return(estimated_abundances)
		else:
			buffer = estimated_abundances
	return(estimated_abundances)


def isoform_abundances(gene_tab, gene):
	'''
	Function to get the isoform relative abundances
	'''
	#check if the tab contain an unique transcript,... so prob = 1
	if len(gene_tab.transcript_name.unique()) == 1:
		return([gene,gene_tab.transcript_name.unique().tolist(),[1]])
	else:
		#go into the em alogorithm
		return [gene,gene_tab.transcript_name.unique().tolist(),em_algorithm(gene_tab).tolist()]


# Fragment file opening
# ---------------------

# print("fragment file opening...")
cell = pd.read_csv(args.cell, sep="\t", names=['bc','gene_name','gene_id','transcript_name','transcript_id','umi','frag_prob_weighted'])

# Perform EM algorithm
print("EM alg...")
res = pd.DataFrame([isoform_abundances(group[1], group[0]) for group in cell[['gene_name','umi','transcript_name','frag_prob_weighted']].groupby("gene_name")])
res.columns = ["gene_name","transcript_name","tr_prob"]
res = res.explode(["transcript_name","tr_prob"])

print("joining...")
cell = cell.drop_duplicates(['gene_name','transcript_name'])
cell = cell.merge(res, on=['gene_name','transcript_name'],  how="left")[['bc','gene_name','gene_id','transcript_name','transcript_id','tr_prob']]

#write into output each cell
print("writing...")
cell.to_csv(args.output_path, sep="\t", index=False)


