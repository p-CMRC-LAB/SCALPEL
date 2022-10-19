
#Author: Franz AKE
#July 2022

#libs
import argparse
import vaex as vx
import pandas as pd
import numpy as np
import dask.dataframe as dd
import sys

#Arguments Parser
parser = argparse.ArgumentParser(description='Processing of Salmon quantification file')
parser.add_argument('-salmon_quant', type=str, default=sys.stdin, help='path of salmon quant file')
parser.add_argument('-output_file', type=str, default=sys.stdout, help='path of output file')
args = parser.parse_args()

#Files
qf = dd.read_csv(args.salmon_quant, sep="\t", header=0, names=['Name','length','EffectiveLength','TPM','NumReads','Sample']).compute()

#expand gene name
name_tr = qf.Name.values[1]
print(name_tr)
if '|' in name_tr:
	qf = pd.concat([qf.Name.str.split("|", expand=True).loc[:,[0,1,4,5]], qf.loc[:,['length','NumReads','Sample']]], axis=1)
	qf.columns = ['transcript_id','gene_id','transcript_name','gene_name','length','NumReads','Sample']
else:
	qf = qf.loc[:,['Name','length','NumReads','Sample']]
	qf.columns = ['transcript_id','length','NumReads','Sample']

#filter transcripts null
qf = qf[~qf.transcript_id.isin(qf[qf.NumReads==0].transcript_id)]
qf = qf.drop_duplicates()
grouped = qf.sort_values(['transcript_id','NumReads']).groupby("transcript_id")
out = list()
for tab in grouped:
	tab[1]['NumReads'] = tab[1]['NumReads'].mean()
	out.append(tab[1])
qf = pd.concat(out)
qf['Sample'] = 'averaged'
qf = qf.drop_duplicates()
#qf = pd.concat([tab[1][tab[1].NumReads==avg(tab[1].NumReads)] for tab in grouped])

#write csv and parquet file
qf.to_csv(args.output_file, sep="\t", doublequote=False, index=False)
