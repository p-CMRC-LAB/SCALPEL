#/usr/bin/python3


#import packages
import argparse
#packages importations
import warnings
warnings.filterwarnings("ignore")
import pandas as pd
import glob
import vaex as vx
import dask.dataframe as dd
import multiprocessing
#
# # Set the Pandas engine to 'numexpr'
# pd.options.mode.chained_assignment = None  # default='warn'
# pd.options.compute.use_numexpr = True



#Argument Parser
#**************
parser = argparse.ArgumentParser(description='Calculate probabilities of fragments...')
parser.add_argument('fragments_path', metavar='Bfip1', type=str, help='path of fragment file')
parser.add_argument('prob_path', metavar='Efip', type=str, help='path of probability bins file')
parser.add_argument('output_bed', metavar='Bmfip', type=str, help='path of output fragment file')
args = parser.parse_args()


#Files opening
#*************
print('File opening...')
reads = vx.open(args.fragments_path)
prob = vx.from_csv(args.prob_path, sep="\t", names=['dist_END', 'counts', 'probs_bin'])


#Merging
#*******
print('File merging...')
reads = reads.join(prob, left_on='dist_END', right_on='dist_END', how='left')
reads['probs_bin'] = reads.probs_bin.fillnan(0)
#reads['probability_normalized'] = reads.probability_normalized.fillnan(0)


#Selection of columns
#********************
print('Column selection...')
reads = reads[['bc','gene_name','gene_id','transcript_name','transcript_id','salmon_rlp','umi','frag_id','probs_bin']]
#scale between 0 and 100
reads['salmon_rlp'] = reads.salmon_rlp * 100


#delete sequencing tags attached to the barcode
reads['bc'] = reads.bc.str.replace("XC:Z:","")
reads['bc'] = reads.bc.str.replace("CB:Z:","")
#delete sequencing tags attached to the umi
reads['umi'] = reads.umi.str.replace("XM:Z:","")
reads['umi'] = reads.umi.str.replace("UB:Z:","")


#Grouping
#********
print('Grouping...')
#reads = reads.astype({'bc':'str','gene_name':'str','gene_id':'str','transcript_name':'str','transcript_id':'str','salmon_rlp':'str','umi':'str','frag_id':'str'})
reads_pd = reads.to_pandas_df()
print(reads)
print("product...")
reads = reads_pd.groupby(['bc','gene_name','gene_id','transcript_name','transcript_id','salmon_rlp','umi','frag_id']).prod().reset_index()


#Splitting table by cells and writing
#************************************
print("writing")
reads.to_csv(args.output_bed, sep="\t", index=False, header=False)
