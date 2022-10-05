
#import packages
import argparse
import pandas as pd
import warnings
import sys
from IPython.display import display
warnings.filterwarnings("ignore")

#Argunent Parser
#**************
parser = argparse.ArgumentParser(description='filtering of intronic exons...')
parser.add_argument('esam', type=str, help='path of exonic sam file')
parser.add_argument('isam', type=str, help='path of intronic sam file')
parser.add_argument('seq_type', type=str, help='sequencing_type')
parser.add_argument('output_file', nargs='?', type=argparse.FileType('w'), default=sys.stdout, help='path of output sam file')
args = parser.parse_args()


#column types
if args.seq_type == 'dropseq':
	coltypes={'read_id':'category','str_rd':'int8','chr_rd':'category','coord':'int64','mq':'category','cigar':'category','c6':'category','c7':'category','c8':'category','read_seq':'str','read_qual':'str','BC':'category','UMI':'category','gfunc_rd':'category'}
elif args.seq_type == 'chromium':
	coltypes={'read_id':'category','str_rd':'int8','chr_rd':'category','coord':'int64','mq':'category','cigar':'category','c6':'category','c7':'category','c8':'category','read_seq':'str','read_qual':'str','BC':'category','UMI':'category'}
else:
	exit('seq_type tag error in bam_filtering.py')

#read intronics and exonics reads
exonics = pd.read_csv(args.esam, sep='\t', header=None)
intronics = {}
intronics['in'] = pd.read_csv(args.isam, sep='\t', header=None)

#identify the columns assciated to barcodes, umi, gene_function_rd, gene_name_rd, strand_mapping and rename
tags_start = 11
for colidx, col in enumerate(exonics.loc[0,tags_start:]):
	if col[:2] in ['XC','CB']:
		exonics = exonics.rename(columns = {tags_start+colidx:'BC'})
	elif col[:2] in ['XM','UB']:
		exonics = exonics.rename(columns = {tags_start+colidx:'UMI'})
	elif col[:2] in ['gf']:
		exonics = exonics.rename(columns = {tags_start+colidx:'gfunc_rd'})
	elif col[:2] in ['gs']:
		exonics = exonics.rename(columns = {tags_start+colidx:'gstrand_rd'})


#names the default bam columns and categorization
exonics = exonics.rename(columns = {0:'read_id', 1:'str_rd', 2:'chr_rd', 3:'coord', 4:'mq', 5:'cigar', 6:'c6', 7:'c7', 8:'c8', 9:'read_seq', 10:'read_qual'})
display(exonics)
exonics = exonics.astype(coltypes)

#name intronics and categorization
intronics['in'].columns = exonics.columns
intronics['in'] = intronics['in'].astype(coltypes)


#process the table if dropseq or chromium
if args.seq_type == 'dropseq':
	#delete sequencing tags in exonics...
	exonics['gfunc_rd'] = exonics['gfunc_rd'].str.replace("gf:Z:","")
	exonics['gstrand_rd'] = exonics['gstrand_rd'].str.replace("gs:Z:","")
	#split and explode...
	exonics['gfunc_rd'] = exonics['gfunc_rd'].str.split(",")
	exonics['gstrand_rd'] = exonics['gstrand_rd'].str.split(",")
	#add index tag
	exonics['index_rd'] = exonics.index
	#explode
	print('exploding...')
	exonics = exonics.explode(['gfunc_rd','gstrand_rd'])
	#select columns
	exonics = exonics[(exonics.gfunc_rd.isin(['CODING','UTR'])) & (exonics.eval('(str_rd == 16 & gstrand_rd == "-") | (str_rd == 0 & gstrand_rd == "+")'))]

	print('duplicating filt...')
	# exonics = pd.concat([tab[1] if len(tab[1]) == 1 else tab[1].loc[:1,:] for tab in exonics.groupby('index_rd')])
	exonics = exonics.drop_duplicates(['index_rd','read_id'])

	#re add sequencing tags
	exonics['gfunc_rd'] = 'gf:Z:' + exonics['gfunc_rd']
	exonics['gstrand_rd'] = 'gs:Z:' + exonics['gstrand_rd']
	del exonics['gfunc_rd']
	del exonics['gstrand_rd']
	del exonics['index_rd']


#filter out intronics associated BC/UMI
exonics['frag_id'] = exonics['BC'].astype(str) + '***' + exonics['UMI'].astype(str)
trs_todisc = intronics['in']['BC'].astype(str) + '***' + intronics['in']['UMI'].astype(str)
intronics.clear()
exonics = exonics[~exonics.frag_id.isin(trs_todisc)]
del exonics['frag_id']

#write file
exonics.to_csv(args.output_file, sep="\t", header=None, index=False)
