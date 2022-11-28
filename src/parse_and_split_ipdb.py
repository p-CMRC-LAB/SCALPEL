


#Script for splitting internal priming reference table according tot the GTF reference annotation (chr or not tagged)
#author: Franz AKE


#libraries importations
import pandas as pd
import argparse
from IPython.display import display


#Argunent Parser
#****************
parser = argparse.ArgumentParser(description='Splitting IPDB reference file')
parser.add_argument('ipdb_file', metavar='IPDB_file_input', type=str, help='path of ipdb file')
parser.add_argument('chr_tag', metavar='chromosome_tag', type=str, help='chromosome tag character')
parser.add_argument('output_ipdb', metavar='IPDB_output', type=str, help='path of filtered ipdb file')
args = parser.parse_args()


#Params
CHR_TAG = args.chr_tag
print(CHR_TAG)

#File opening
print("file opening..")
ipdb = pd.read_csv(args.ipdb_file, sep="\t", header=None, names = ['chr_ipn','start_ipn','end_ipn','c3','c4','str_ipn'], dtype = {'chr_ipn':'str'})

#filter ipdb table
ipdb = ipdb[['chr_ipn','start_ipn','end_ipn','str_ipn']]
ipdb.chr_ipn = ipdb.chr_ipn.str.replace("chr","")

if len(CHR_TAG) > 3:
	#assume that we have the presence of a chr character in the chromosome col --> remove it
	CHR_TAG = CHR_TAG[3:]
else:
	pass

#filtering IPDB table
ipdb = ipdb[ipdb.chr_ipn == CHR_TAG].sort_values(['chr_ipn','start_ipn','end_ipn','str_ipn']).reset_index(drop=True)

#write table
ipdb.to_csv(args.output_ipdb, sep="\t", header=False, index=False)
