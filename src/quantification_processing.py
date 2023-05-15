
#Author: Franz AKE
#July 2022

#libs
import argparse
import pandas as pd

#Arguments Parser
parser = argparse.ArgumentParser(description='Processing of Salmon quantification file')
parser.add_argument('-salmon_quant', type=str, help='path of salmon quant file')
parser.add_argument('-output_file', type=str, help='path of output file')
args = parser.parse_args()

#Files
print("read file...")
all_samples = []
for filin in args.salmon_quant.split(","):
	#open file
	filetab = pd.read_csv(filin, sep="\t")
	filetab["Sample"] = filin
	all_samples.append(filetab)
qf = pd.concat(all_samples)

#expand gene name
name_tr = qf.Name.values[1]
if '|' in name_tr:
	qf = pd.concat([qf.Name.str.split("|", expand=True).loc[:,[0]], qf[['Length', 'NumReads','Sample']]], axis=1)
else:
	qf = qf.loc[:,['Name','Length','NumReads','Sample']]
qf.columns = ['transcript_id','length','NumReads','Sample']


#filter transcripts null
qf = qf.drop_duplicates()
grouped = qf.sort_values(['transcript_id','NumReads']).groupby("transcript_id")
out = list()
for tab in grouped:
	tab[1]['NumReads'] = tab[1]['NumReads'].mean()
	out.append(tab[1])
qf = pd.concat(out)
qf = qf[~qf.transcript_id.isin(qf[qf.NumReads==0].transcript_id)]
qf['Sample'] = 'averaged'
qf = qf.drop_duplicates()

#write csv and parquet file
qf.to_csv(args.output_file, sep="\t", doublequote=False, index=False)
