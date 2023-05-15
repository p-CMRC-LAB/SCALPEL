


from gtfparse import read_gtf
import pandas as pd
import argparse


#Argument Parser
#****************
parser = argparse.ArgumentParser(description='Processing of GTF file')
parser.add_argument('gtf_file', metavar='GTF_file_input', type=str, help='path of gtf file')
parser.add_argument('quant_file', metavar='quant_file_input', type=str, help='path of transcript quantification file by salmon')
args = parser.parse_args()



#Files reading
#*************
print('file opening')
gtf = read_gtf(args.gtf_file)
qf = pd.read_csv(args.quant_file, sep="\t").dropna()



# (1) Check biotype/type concordance in gtf colnames (gencode annotation differences...)
if 'gene_biotype' in gtf.columns:
	gtf = gtf.rename(columns = {'gene_biotype':'gene_type', 'transcript_biotype':'transcript_type'})
# (2) Check transcript_version and gene_id version (UCSC name)
if ('gene_version' in gtf.columns) and ('transcript_version' in gtf.columns):
	gtf['gene_id'] = gtf.gene_id + '.' + gtf.gene_version
	gtf['transcript_id'] = gtf.transcript_id + '.' + gtf.transcript_version

# (3) Filter transcripts to specific gene type...
gtf = gtf[(gtf.feature.isin(['exon','UTR'])) & (gtf.gene_type.isin(['protein_coding','processed_pseudogene','transcribed_unprocessed_pseudogene','lncRNA']))]
gtf = gtf[["seqname","feature","start","end","strand","gene_id","gene_name","gene_type","transcript_id","transcript_name","transcript_type","exon_number","exon_id"]]


# (4) Filter gtf file base on Salmon quantification
print("Filter gtf file based on Salmon quantification...")
# add gene_name annotation in Salmon quant
qf = pd.merge(qf, gtf[['transcript_id','gene_name']], on='transcript_id', how='left').drop_duplicates().dropna()
gtf = gtf[gtf.transcript_id.isin(qf.transcript_id)]

#function to calculate average quantification in isoform/gene
def rlp_transcripts(tab):
	tab['salmon_rlp'] = tab.NumReads / tab.NumReads.sum()
	return(tab)
qf = qf.groupby('gene_name').apply(lambda x: rlp_transcripts(x))

#filtering
gtf = gtf[gtf.transcript_id.isin(qf.transcript_id)]
gtf = pd.merge(gtf, qf[['transcript_id','length','NumReads','salmon_rlp']], on=['transcript_id'], how='left').drop_duplicates().dropna()

# Also add +1 to match the index format as pyrange bedify the GTF input
# gtf.Start = gtf.Start + 1


# (5) write chromosome files
print("writing...")
seqnames = gtf.seqname.unique()
for chr in seqnames:
	gtf[gtf.seqname==chr].to_csv(chr, sep="\t", index=False)

