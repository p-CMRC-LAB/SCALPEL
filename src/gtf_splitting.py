

#library importation
#*******************
import pyranges as pr
import pandas as pd
import argparse
import numpy as np
from IPython.display import display


#Argunent Parser
#****************
parser = argparse.ArgumentParser(description='Processing of GTF file')
parser.add_argument('gtf_file', metavar='GTF_file_input', type=str, help='path of gtf file')
parser.add_argument('quant_file', metavar='quant_file_input', type=str, help='path of transcript quantification file by salmon')
args = parser.parse_args()


#Files reading
#*************
print('file opening')
gtf = pr.read_gtf(args.gtf_file, as_df=True)
qf = pd.read_csv(args.quant_file, sep="\t")

#Main
#****

# (1) Check biotype/type concordance in gtf colnames (gencode annotation differences...)
if 'gene_biotype' in gtf.columns:
	gtf = gtf.rename(columns = {'gene_biotype':'gene_type', 'transcript_biotype':'transcript_type'})
else:
	pass

# (2) Check transcript_version and gene_id version (UCSC name)
if ('gene_version' in gtf.columns) and ('transcript_version' in gtf.columns):
	gtf['gene_id'] = gtf.gene_id + '.' + gtf.gene_version
	gtf['transcript_id'] = gtf.transcript_id + '.' + gtf.transcript_version

# (3) Filter transcripts to specific gene type...
print('gtf opening...')
gtf = gtf[(gtf.Feature=='exon') & (gtf.gene_type.isin(['protein_coding','processed_pseudogene','transcribed_unprocessed_pseudogene','lncRNA']))]
gtf = gtf[['Chromosome','Start','End','Strand','gene_id','gene_name','gene_type','transcript_id','transcript_name','transcript_type','exon_number','exon_id']]

# (4) Compute transcripts salmon relative probabilities...
if 'gene_name' not in  qf.columns:
	qf = pd.merge(qf, gtf[['transcript_id','gene_name']], on='transcript_id', how='right').drop_duplicates()
	qf = qf.dropna()

def rlp_transcripts(tab):
	tab['salmon_rlp'] = tab.NumReads / tab.NumReads.sum()
	return(tab)
qf = qf.groupby('gene_name').apply(lambda x: rlp_transcripts(x))

# (5) Merge salmon quantification...
gtf = pd.merge(gtf, qf[['transcript_id','length','NumReads','salmon_rlp']], on=['transcript_id'], how='right').drop_duplicates()

# (6) Discard Ribosomical genes...
genes_names = gtf['gene_name']
# genes_todelete = (genes_names[genes_names.str.contains("^RP[L|S]")]).tolist()
# gtf = (gtf[~gtf.gene_name.isin(genes_todelete)])

# (7) Discard Mitochondrial chromosome associated...
gtf = gtf[(gtf.Chromosome != 'chrM') & (gtf.Chromosome != 'M')]
gtf = gtf[(gtf.Chromosome != 'chrMT') & (gtf.Chromosome != 'MT')]

# (8) Get chromosome names list...
gtf = gtf.dropna()
chromosome_names = gtf.Chromosome.unique()

# Also add +1 to match the index format as pyrange bedify the GTF input
gtf.Start = gtf.Start + 1

# (9) Order gtf file...
gtf = gtf.sort_values(by=['Chromosome','Start','End'], ascending=[True,True,False]).reset_index(drop=True)

# (10) Split transcript GTF file...
print('writing...')
for chr in chromosome_names:
	gtf[gtf.Chromosome==chr].to_csv(chr, sep="\t", index=False)
