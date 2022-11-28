
#import packages
import subprocess
import csv
import argparse
import warnings
warnings.filterwarnings("ignore")
import pandas as pd
import gc
from IPython.display import display

#Argunent Parser
#**************
parser = argparse.ArgumentParser(description='Overlapping of reads with Exons')
parser.add_argument('bedf', metavar='reads', type=str, help='path of reads file1')
parser.add_argument('exons', metavar='exons', type=str, help='path of exons file1')
parser.add_argument('bed_gtf_path', metavar='Bfip1', type=str, help='path of bed (bedtools) file1')
parser.add_argument('output_path', metavar='Bmfip', type=str, help='path of merge bed_gtf file')
args = parser.parse_args()

#Params
#******
MARGIN_END = 10

#files opening
print('file opening...')
#read files
bed_gtfs = pd.read_csv(args.bed_gtf_path, sep="\*\*\*\*", header=None, names=['reads','exons'], dtype={'reads': 'object', 'exons': 'category'})
beds_all = pd.read_parquet(args.bedf)
beds_all['rd_id'] = beds_all.index.astype('int32')
beds_all['frag_id'] = beds_all.bc.astype('str') + '/' + beds_all.umi.astype('str')

#get the number of frag_id
container = {}
container['in'] = beds_all.groupby('frag_id').size().reset_index()
container['in'].columns = ['frag_id','frag_id_nb']
beds_all = beds_all.merge(container['in'], on='frag_id', how='left')
container.clear()

#exons file
exons = pd.read_parquet(args.exons)
exons['ex_id'] = exons.index.astype('int32')

#get unmapped and mapped reads and exons
unmappeds = bed_gtfs[bed_gtfs.exons.isnull()]
mappeds = bed_gtfs[~bed_gtfs.exons.isnull()]

display(bed_gtfs)
display(mappeds)
display(unmappeds)

if len(mappeds) != 0:
	#let's split the diferent exons mappeds
	mappeds['exons'] = pd.DataFrame(mappeds.exons.str.split(';'))
	mappeds = mappeds.explode('exons')
	print('writing...')
	mappeds.to_csv('mappeds', sep='\t', doublequote=False, index=False, header=False, quoting = csv.QUOTE_NONE, escapechar = ',')
	subprocess.run(' sed "s/\,//g" mappeds > mappeds2', shell=True)
	print('reading...')
	mappeds = pd.read_csv('mappeds2', sep='\t', names=['chr_rd','start_rd','end_rd','rd_id','Chromosome','Start','End','ex_id'], dtype={'Chromosome':'str', 'chr_rd':'category','ex_id':'int32'})

	#let's add exons information
	mappeds = mappeds.merge(exons, on=['Chromosome','Start','End','ex_id'], how='left')
	#read_information
	mappeds = mappeds.merge(beds_all, on=['chr_rd','start_rd','end_rd','rd_id'], how='left')

	#let's check bordering conditions
	print('filtering border conditions...')
	#checkings
	mappeds = mappeds[mappeds.str_rd==mappeds.Strand][['chr_rd','start_rd','end_rd','str_rd','Start','End','StartR','EndR','frag_id','read_id','bc','umi','gene_id','gene_name','transcript_id','transcript_name','exon_id','exon_number','nb_splices','salmon_rlp','frag_id_nb']]
	#discard fragment overlapping out of authorized area
	mappeds['frag_id_transcript'] = mappeds.frag_id.astype('str') + "-" + mappeds.transcript_name.astype('str')
	#coding areas
	to_discard = mappeds[(mappeds.exon_number!=1) & ((mappeds.start_rd < mappeds.Start) | (mappeds.end_rd > mappeds.End))][['frag_id','transcript_name']]
	to_discard = to_discard.frag_id.astype('str') + '-' + to_discard.transcript_name.astype('str')
	mappeds = mappeds[~mappeds.frag_id_transcript.isin(to_discard)]
	#utr areas
	to_discard = mappeds[(mappeds.exon_number==1) & (((mappeds.str_rd=='+') & (mappeds.start_rd < mappeds.Start)) | ((mappeds.str_rd=='-') & (mappeds.end_rd > mappeds.End)))][['frag_id','transcript_name']]
	to_discard = to_discard.frag_id.astype('str') + '-' + to_discard.transcript_name.astype('str')
	mappeds = mappeds[~mappeds.frag_id_transcript.isin(to_discard)]
	#utrs ends
	to_discard = mappeds[(mappeds.exon_number==1) & (((mappeds.str_rd=='+') & (mappeds.end_rd > mappeds.End + MARGIN_END)) | ((mappeds.str_rd=='-') & (mappeds.start_rd < mappeds.Start - MARGIN_END)))][['frag_id','transcript_name']]
	to_discard = to_discard.frag_id.astype('str') + '-' + to_discard.transcript_name.astype('str')
	mappeds = mappeds[~mappeds.frag_id_transcript.isin(to_discard)]

	#drop duplicates
	mappeds = mappeds.drop_duplicates()

	#filter fragments associated to unmappeds reads
	unmappeds = unmappeds.reads.str.split('\t', expand=True)
	unmappeds.columns = ['chr_rd','start_rd','end_rd','rd_id']
	unmappeds = unmappeds.astype({'chr_rd':'str','start_rd':'int64','end_rd':'int64','rd_id':'int64'})
	unmappeds = unmappeds.merge(beds_all, on=['chr_rd','start_rd','end_rd','rd_id'], how='left')
	unmappeds = unmappeds.bc.astype('str') + '/' + unmappeds.umi.astype('str')
	mappeds = mappeds[~mappeds.frag_id.isin(unmappeds)]

	print('filtering frag_id number...')
	#filter based on unaccordance in fragment number
	container = {}
	container['in'] = mappeds.groupby(['frag_id_transcript','frag_id_nb']).size().reset_index()
	container['in'].columns = ['frag_id_transcript','frag_id_nb','frag_id_nb_obs']
	container['in'] = container['in'][container['in'].frag_id_nb==container['in'].frag_id_nb_obs]
	mappeds = mappeds[mappeds.frag_id_transcript.isin(container['in'].frag_id_transcript)].reset_index(drop=True)
	del mappeds['frag_id_transcript']

	#Writing
	mappeds = mappeds.astype({'frag_id':'category', 'frag_id_nb':'int8', 'read_id':'category', 'bc':'category', 'umi':'category'})
	mappeds.to_parquet(args.output_path, engine='pyarrow')

else:
	print("Empty mapped reads on the transcriptome reference... Check about it!")



# #Params
# #******
# coltypes = {'chr_rd':'category','start_rd':'int64','end_rd':'int64','str_rd':'category','read_id':'category','bc':'category','umi':'category','read_id':'str','nb_splices':'int8'}
# ex_coltypes = {'Chromosome':'str','Start':'int64','End':'int64','Strand':'str','StartR':'int32','EndR':'int32','gene_id':'category','gene_name':'category','transcript_id':'category','transcript_name':'category','exon_id':'category','exon_number':'int8','gene_type':'category','transcript_type':'category','salmon_rlp':'float'}
# MARGIN_END = 7
#
#
# #Files opening
# print("Files opening...")
# beds_all = {}
# bed_gtfs = {}
# bed_gtfs['in'] = pd.read_csv(args.bed_gtf_path, sep="\*\*\*\*", header=None, names=['reads','exons'], dtype={'reads': 'object', 'exons': 'category'})
# beds_all['in'] = pd.read_csv(args.bedf, sep="\t", header=None, names = ['chr_rd','start_rd','end_rd','read_id','str_rd','bc','umi','nb_splices'], dtype = coltypes)
# exons = pd.read_csv(args.exons, sep="\t", dtype = ex_coltypes)
#
# #get unmapped and mapped reads
# print("get unmapped and mapped reads...")
# unmappeds = {}
# unmappeds['in'] = bed_gtfs['in'][bed_gtfs['in'].exons.isnull()]
# bed_gtfs['in'] = bed_gtfs['in'][~bed_gtfs['in'].exons.isnull()]
# #get unmappeds['in'] frag_id infos
# if len(unmappeds['in']) != 0:
# 	unmappeds['in'] = unmappeds['in'].reads.str.split('\t', expand = True)
# 	unmappeds['in'].columns = ['chr_rd','start_rd','end_rd','read_id']
# 	unmappeds['in'] = unmappeds['in'].astype({'chr_rd':'category','start_rd':'int64','end_rd':'int64','read_id':'str'})
# 	unmappeds['in'] = unmappeds['in'].merge(beds_all['in'][['chr_rd','start_rd','end_rd','read_id','bc','umi']], how='left')
# 	unmappeds['in'] = unmappeds['in'].bc.astype(str) + '/' + unmappeds['in'].umi.astype(str)
#
# #process mapped reads
# print("get unmapped and mapped reads...")
# if len(bed_gtfs['in']) != 0:
# 	#reads
# 	print("reads...")
# 	mappeds = bed_gtfs['in'].reads.str.split('\t', expand=True)
# 	mappeds.columns = ['chr_rd','start_rd','end_rd','read_id']
# 	mappeds = mappeds.astype({'chr_rd':'category','start_rd':'int64','end_rd':'int64','read_id':'str'})
# 	mappeds = mappeds.astype({'chr_rd':'category','start_rd':'int64','end_rd':'int64','read_id':'str'})
# 	mappeds['idx_map'] = mappeds.index
# 	mappeds = mappeds.merge(beds_all['in'][['chr_rd','start_rd','end_rd','read_id','str_rd','bc','umi','nb_splices']], how='left')
# 	mappeds['frag_id'] = mappeds.bc.astype(str) + '/' + mappeds.umi.astype(str)
# 	beds_all.clear()
#
# 	#discard all the bc/umi associated to umappeds
# 	mappeds = mappeds[~mappeds.frag_id.isin(unmappeds['in'])]
# 	unmappeds.clear()
#
# 	#calculate frag_id number
# 	print("frag_id number...")
# 	container = {}
# 	container['in'] = mappeds.groupby('frag_id').size().reset_index()
# 	container['in'].columns = ['frag_id','frag_id_nb']
# 	mappeds = mappeds.merge(container['in'], on='frag_id', how='left')
# 	container.clear()
# 	display(mappeds.info(memory_usage=True))
#
# 	#exons
# 	print("exons...")
# 	exons_mappeds = {}
# 	exons_mappeds['in'] = pd.DataFrame(bed_gtfs['in'].exons.str.split(';'))
# 	bed_gtfs.clear()
#
# 	exons_mappeds['in']['idx'] = exons_mappeds['in'].index
# 	print('exploding...')
# 	exons_mappeds['in'] = (exons_mappeds['in'].explode('exons'))
# 	display(exons_mappeds['in'].info(memory_usage=True))
# 	# exons_mappeds['in'] = exons_mappeds['in'].drop_duplicates(['exons','idx'])
# 	edix = exons_mappeds['in'].idx.tolist()
#
# 	#let's write this exon first split in a file and read it again with the right casting of variable
#
# 	display('expanding...')
# 	exons_mappeds['in'] = exons_mappeds['in'].exons.str.split('\t', expand=True)
# 	exons_mappeds['in'].columns = ['Chromosome','Start','End','gene_name','transcript_name','exon_number']
# 	exons_mappeds['in']['idx_map'] = edix
# 	exons_mappeds['in'] = exons_mappeds['in'].astype({'Chromosome':'str','Start':'int64','End':'int64','gene_name':'category','transcript_name':'category','exon_number':'int8'})
#
# 	print('merging...')
# 	exons_mappeds['in'] = exons_mappeds['in'].merge(exons, how='left', on = ['Chromosome','Start','End','gene_name','transcript_name','exon_number'])
#
# 	#join tables
# 	print("join tables...")
# 	mappeds = mappeds.merge(exons_mappeds['in'], on='idx_map', how='left')
# 	exons_mappeds.clear()
#
# 	#checkings
# 	mappeds = mappeds[mappeds.str_rd==mappeds.Strand][['chr_rd','start_rd','end_rd','str_rd','Start','End','StartR','EndR','frag_id','read_id','bc','umi','gene_id','gene_name','transcript_id','transcript_name','exon_id','exon_number','nb_splices','salmon_rlp','frag_id_nb']]
# 	#discard fragment overlapping out of authorized area
# 	mappeds['frag_id_transcript'] = mappeds.frag_id.astype('str') + "-" + mappeds.transcript_name.astype('str')
# 	#coding areas
# 	to_discard = mappeds[(mappeds.exon_number!=1) & ((mappeds.start_rd < mappeds.Start) | (mappeds.end_rd > mappeds.End))][['frag_id','transcript_name']]
# 	to_discard = to_discard.frag_id + '-' + to_discard.transcript_name
# 	mappeds = mappeds[~mappeds.frag_id_transcript.isin(to_discard)]
# 	#utr areas
# 	to_discard = mappeds[(mappeds.exon_number==1) & (((mappeds.str_rd=='+') & (mappeds.start_rd < mappeds.Start)) | ((mappeds.str_rd=='-') & (mappeds.end_rd > mappeds.End)))][['frag_id','transcript_name']]
# 	to_discard = to_discard.frag_id + '-' + to_discard.transcript_name
# 	mappeds = mappeds[~mappeds.frag_id_transcript.isin(to_discard)]
# 	#utrs ends
# 	to_discard = mappeds[(mappeds.exon_number==1) & (((mappeds.str_rd=='+') & (mappeds.end_rd > mappeds.End + MARGIN_END)) | ((mappeds.str_rd=='-') & (mappeds.start_rd < mappeds.Start - MARGIN_END)))][['frag_id','transcript_name']]
# 	to_discard = to_discard.frag_id + '-' + to_discard.transcript_name
# 	mappeds = mappeds[~mappeds.frag_id_transcript.isin(to_discard)]
#
#
# 	#frag_id_nb_filtering
# 	print('checking...')
# 	mappeds = mappeds.astype({'frag_id':'str','transcript_name':'str','frag_id_nb':'int8'})
#
# 	container = {}
# 	container['in'] = mappeds.groupby(['frag_id','transcript_name','frag_id_nb']).size().reset_index()
# 	container['in'].columns = ['frag_id','transcript_name','frag_id_nb','frag_id_nb_obs']
# 	container['in'] = container['in'][container['in'].frag_id_nb==container['in'].frag_id_nb_obs]
# 	container['in']['infilter'] = container['in'].frag_id + '-' + container['in'].transcript_name + '-' + container['in'].frag_id_nb.astype(str)
# 	mappeds['infilter'] = mappeds.frag_id + '-' + mappeds.transcript_name + '-' + mappeds.frag_id_nb.astype(str)
# 	mappeds = mappeds[mappeds.infilter.isin(container['in'].infilter)]
# 	container.clear()
# 	del mappeds['infilter']
# 	del mappeds['frag_id_transcript']
#
# 	#Writing
# 	mappeds = mappeds.astype({'gene_name':'category', 'frag_id_nb':'int8'})
# 	mappeds.to_parquet(args.output_path, engine='pyarrow')
#
# else:
# 	exit("Empty mapped reads on the transcriptome reference... Check about it!")
