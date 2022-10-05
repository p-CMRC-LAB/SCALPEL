

#import packages
import argparse
import warnings
warnings.filterwarnings("ignore")
import pandas as pd
import gc
import numpy as np
from IPython.display import display
import dask.dataframe as dd


#Argument Parser
#**************
parser = argparse.ArgumentParser(description='Filtering of SAM file')
parser.add_argument('sam_path', metavar='Bfip1', type=str, help='path of sam file')
parser.add_argument('fragments_path', metavar='Bfip2', type=str, help='path of fragments file')
parser.add_argument('output_sam_path', metavar='Bmfip', type=str, help='path of filtered sam file')
args = parser.parse_args()


#open SAM file & Fragment files
print('File opening...')
samf = dd.read_csv(args.sam_path, sep="\t", header=None).compute()
samf = samf.rename(columns = {0:'read_id', 1:'str_rd', 2:'chr_rd', 3:'coord', 4:'mq', 5:'cigar', 6:'c6', 7:'c7', 8:'c8', 9:'read_seq', 10:'read_qual', 11:'bc', 12:'umi'})

fragments = dd.read_parquet(args.fragments_path, engine='pyarrow').compute().drop_duplicates()
fragments.bc = 'XC:Z:' + fragments.bc.astype(str)
fragments.umi = 'XM:Z:' + fragments.umi.astype(str)

#filtering
samf = samf[(samf.read_id.isin(fragments.read_id)) & (samf.bc.isin(fragments.bc)) & (samf.umi.isin(fragments.umi))]

#Writing
print('Writing...')
samf.to_csv(args.output_sam_path, sep="\t", header=False, index=False, doublequote=False)




#
