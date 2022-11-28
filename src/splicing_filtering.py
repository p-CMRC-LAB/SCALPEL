

#import packages
import argparse
import pandas as pd
import warnings
import sys
from IPython.display import display
warnings.filterwarnings("ignore")

#Argunent Parser
#***************
parser = argparse.ArgumentParser(description='preprocess spliced reads...')
parser.add_argument('bed', nargs='?', type=argparse.FileType('r'), default=sys.stdin, help='path of bed file')
parser.add_argument('output_file', type=str, help='path of output bed file')
parser.add_argument('bedmap_output', nargs='?', type=argparse.FileType('w'), default=sys.stdout, help='path of output bedmap format file')
args = parser.parse_args()



#Files opening
coltypes = {'chr_rd':'category', 'start_rd':'int64','end_rd':'int64','read_id':'str','str_rd':'category','bc':'category','umi':'category'}
reads = pd.read_csv(args.bed, sep='\t', names = ['chr_rd','start_rd','end_rd','read_id','mq','str_rd','c6','c7','c8','c9','c10','c11','c12','bc','umi'], dtype = coltypes, usecols = ['chr_rd','start_rd','end_rd','read_id','str_rd','bc','umi'])


#delete seq tags of barcodes, umis
#dropseq
reads.bc = reads.bc.str.replace('XC:Z:','')
reads.umi = reads.umi.str.replace('XM:Z:','')
#chromium
reads.bc = reads.bc.str.replace('CB:Z:','')
reads.umi = reads.umi.str.replace('UB:Z:','')

#tag reads with effective splicing
reads['read_id2'] = reads.read_id.str.replace("/[0-9]","")

container = {}
container['in'] = pd.DataFrame([[tab[0],len(tab[1])] for tab in reads.groupby('read_id2')])
container['in'].columns = ['read_id2','nb_splices']
reads = reads.merge(container['in'], on='read_id2', how='left')
container.clear()

reads = reads.sort_values(['chr_rd','start_rd','end_rd']).reset_index(drop=True)
reads.read_id = reads.read_id2
del reads['read_id2']

#add +1 to the start coordinates to match with the gtf file coordinates
# reads.start_rd = reads.start_rd + 1

# #writing
reads.to_parquet(args.output_file, engine='pyarrow')
reads['id'] = reads.index
reads[['chr_rd','start_rd','end_rd','id']].to_csv(args.bedmap_output, sep="\t", header=None, index=False)
