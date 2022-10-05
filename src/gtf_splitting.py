

#libs
import os
import vaex as vx
import pyranges as pr
import pandas as pd
import argparse
import numpy as np
from IPython.display import display



#Argunent Parser
#**************
parser = argparse.ArgumentParser(description='Processing of GTF file')
parser.add_argument('gtf_file', metavar='GTF_file_input', type=str, help='path of gtf file splitted by chromosome')
parser.add_argument('quant_file', metavar='quant_file_input', type=str, help='path of transcript quantification file by salmon')
args = parser.parse_args()



QF_FLAG = True





gtf = pr.read_gtf(args.gtf_file, as_df=True)

if args.quant_file != 'null':
	qf = pd.read_csv(args.quant_file, sep="\t")
else:
	QF_FLAG = False

#check biotype/type concordance in gtf colnames
if 'gene_biotype' in gtf.columns:
	gtf = gtf.rename(columns = {'gene_biotype':'gene_type', 'transcript_biotype':'transcript_type'})
else:
	pass

#check transcript_version and gene_id version
if ('gene_version' in gtf.columns) and ('transcript_version' in gtf.columns):
	gtf['gene_id'] = gtf.gene_id + '.' + gtf.gene_version
	gtf['transcript_id'] = gtf.transcript_id + '.' + gtf.transcript_version

# check presence of chr in gtf
if 'chr' in gtf.Chromosome.values[1]:
	gtf['Chromosome'] = gtf.Chromosome.str.replace('chr','')
	# pass
else:
	pass

#extract exons
print('gtf opening...')
gtf = gtf[(gtf.Feature=='exon') & (gtf.gene_type.isin(['protein_coding','processed_pseudogene','transcribed_unprocessed_pseudogene','lncRNA']))]
display(gtf.columns)
gtf = gtf[['Chromosome','Start','End','Strand','gene_id','gene_name','gene_type','transcript_id','transcript_name','transcript_type','exon_number','exon_id']]

display(gtf)
#compute transcripts salmon relative probabilities
def rlp_transcripts(tab):
	tab['salmon_rlp'] = tab.NumReads / tab.NumReads.sum()
	return(tab)
qf = qf.groupby('gene_name').apply(lambda x: rlp_transcripts(x))

#merge salmon quantification
if QF_FLAG == True:
	gtf = pd.merge(gtf, qf[['transcript_name','length','NumReads','salmon_rlp']], on=['transcript_name'], how='inner')
else:
	gtf['length'] = np.array([-1] * len(gtf))
	gtf['NumReads'] = np.array([-1] * len(gtf))

display(gtf)

#Discard Ribosomical genes
genes_names = gtf['gene_name']
genes_todelete = (genes_names[genes_names.str.contains("^RP[L|S]")]).tolist()
gtf = (gtf[~gtf.gene_name.isin(genes_todelete)])

#Discard Mtichondrial chromosome
gtf = gtf[(gtf.Chromosome != 'chrM') & (gtf.Chromosome != 'M')]
gtf = gtf[(gtf.Chromosome != 'chrMT') & (gtf.Chromosome != 'MT')]

#optional (remove chr)
# gtf.Chromosome = gtf.Chromosome.str.replace('chr','')

#get chromosome names list
chromosome_names = gtf.Chromosome.unique()

# Order gtf file
gtf = gtf.sort_values(by=['Chromosome','Start','End'], ascending=[True,True,False]).reset_index(drop=True)

#loop on table and write
print('writing...')
for chr in chromosome_names:
	gtf[gtf.Chromosome==chr].to_csv(chr, sep="\t", index=False)
