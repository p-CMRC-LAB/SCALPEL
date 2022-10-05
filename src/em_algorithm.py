

import pandas as pd
import numpy as np
import argparse
from IPython.display import display


#Argument Parser
#**************
parser = argparse.ArgumentParser(description='compute EM algorithm...')
parser.add_argument('cell_path', metavar='Bfip1', type=str, help='path of cell file')
parser.add_argument('output_path', metavar='Efip', type=str, help='path of probability transcripts file')
args = parser.parse_args()


#Functions
#*********
def posterior(probs, tr_estimates):
	#function to impute transcript relative estimate abundances
	numerator = (probs * tr_estimates)
	#display(numerator)
	#display(sum(numerator))
	#return(numerator/sum(numerator))
	return(numerator.div(numerator.sum(axis=1).values[0]))

def em_algorithm(tab, max_iteration=100):
	#function to comupte em algorithm on transcript table
	if len(tab) == 0:
		return np.nan
	else:
		#get initial estimated abundances starter
		estimated_abundances = [1/tab.shape[1]] * tab.shape[1]
		buffer = estimated_abundances
		#Iterate
		for i in range(max_iteration):
			#get posterior probs
			#display(tab)
			#post_probs = tab.apply(lambda row: posterior(row, estimated_abundances), axis=1)
			post_probs = tab.groupby('umi').apply(lambda x: posterior(x, estimated_abundances))
			#get transcript relative coord estimate
			estimated_abundances = np.round((post_probs.mean(axis=0).values),2)
			#keep track if different else stop
			if False not in (buffer == estimated_abundances):
				return(estimated_abundances.tolist())
			else:
				buffer = estimated_abundances
		return(estimated_abundances.tolist())


#Cell file opening
cell = pd.read_csv(args.cell_path, sep='\t')
cell['frag_prob_weighted'] = cell.probs_bin * cell.salmon_rlp.astype('float') * 1e5
cell = cell.astype({'gene_name':'str'})

#EM algorithm
print('looping')
res = [[tab[1].transcript_name.unique().tolist(), em_algorithm(tab[1][['umi','transcript_name','frag_prob_weighted']].pivot_table(index='umi',columns='transcript_name',values='frag_prob_weighted', fill_value=0))] for tab in cell.groupby('gene_name')]
res = pd.DataFrame(res).dropna()
display(res)
res.columns = ['transcript_name','tr_prob']
res = res.explode(['transcript_name','tr_prob'])
res = cell.merge(res, on='transcript_name', how='left')

#writing
res.to_csv(args.output_path, sep='\t', index=False, doublequote=False)
