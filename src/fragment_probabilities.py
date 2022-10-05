#/usr/bin/python3


#import packages
import dask.dataframe as dd
import argparse
#packages importations
import warnings
warnings.filterwarnings("ignore")
import pandas as pd
import glob



#Argument Parser
#**************
parser = argparse.ArgumentParser(description='Calculate probabilities of fragments...')
parser.add_argument('fragments_path', metavar='Bfip1', type=str, help='path of fragment file')
parser.add_argument('prob_path', metavar='Efip', type=str, help='path of probability bins file')
args = parser.parse_args()


#Files opening
#*************
print('File opening...')
all_files = args.fragments_path
print(all_files)
reads = pd.read_parquet(all_files)
print(reads)
prob = pd.read_csv(args.prob_path, sep="\t")
print(prob)
prob.columns = ['dist_END', 'counts', 'probs_bin', 'counts_normalized', 'probability_normalized']

#Merging
#*******
print('File merging...')
reads = reads.merge(prob, on='dist_END', how='left')
reads[['probs_bin']] = reads[['probs_bin']].fillna(value=0)
reads[['probability_normalized']] = reads[['probability_normalized']].fillna(value=0)

#Selection of columns
#********************
print('Column selection...')
reads = reads[['bc','gene_name','gene_id','transcript_name','transcript_id','salmon_rlp','umi','frag_id','probs_bin']]
#scale between 0 and 100
reads.salmon_rlp = reads.salmon_rlp * 100

#Grouping
#********
print('Grouping...')
reads = reads.astype({'bc':'str','gene_name':'str','gene_id':'str','transcript_name':'str','transcript_id':'str','salmon_rlp':'str','umi':'str','frag_id':'str'})
reads = reads.groupby(['bc','gene_name','gene_id','transcript_name','transcript_id','salmon_rlp','umi','frag_id']).prod().reset_index()

#ponderate fragments probability with salmon weights
#reads.probs_bin = reads.salmon_rlp.astype(float) * reads.probs_bin


#Splitting table by cells and writing
#************************************
print('Splitting and writing...')
print(reads)
for cell in reads.groupby('bc'):
	cell[1].to_csv(cell[0] + '.cell', sep="\t", index=False)
