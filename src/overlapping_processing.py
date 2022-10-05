
#import packages
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
coltypes = {'chr_rd':'category','start_rd':'int64','end_rd':'int64','str_rd':'category','read_id':'category','bc':'category','umi':'category','read_id':'str','nb_splices':'int8'}
ex_coltypes = {'Chromosome':'str','Start':'int64','End':'int64','Strand':'str','StartR':'int32','EndR':'int32','gene_id':'category','gene_name':'category','transcript_id':'category','transcript_name':'category','exon_id':'category','exon_number':'int8','gene_type':'category','transcript_type':'category','salmon_rlp':'float'}
MARGIN_END = 7


#Files opening
bed_gtfs = pd.read_csv(args.bed_gtf_path, sep="\*\*\*\*", header=None, names=['reads','exons'], dtype={'reads': 'object', 'exons': 'category'})
beds_all = pd.read_csv(args.bedf, sep="\t", header=None, names = ['chr_rd','start_rd','end_rd','read_id','str_rd','bc','umi','nb_splices'], dtype = coltypes)
exons = pd.read_csv(args.exons, sep="\t", dtype = ex_coltypes)

#get unmapped and mapped reads
unmappeds = bed_gtfs[bed_gtfs.exons.isnull()]
bed_gtfs = bed_gtfs[~bed_gtfs.exons.isnull()]
#get unmappeds frag_id infos
if len(unmappeds) != 0:
	unmappeds = unmappeds.reads.str.split('\t', expand = True)
	unmappeds.columns = ['chr_rd','start_rd','end_rd','read_id']
	unmappeds = unmappeds.astype({'chr_rd':'category','start_rd':'int64','end_rd':'int64','read_id':'str'})
	unmappeds = unmappeds.merge(beds_all[['chr_rd','start_rd','end_rd','read_id','bc','umi']], how='left')
	unmappeds = unmappeds.bc.astype(str) + '/' + unmappeds.umi.astype(str)

#process mapped reads
if len(bed_gtfs) != 0:
	#reads
	mappeds = bed_gtfs.reads.str.split('\t', expand=True)
	mappeds.columns = ['chr_rd','start_rd','end_rd','read_id']
	mappeds = mappeds.astype({'chr_rd':'category','start_rd':'int64','end_rd':'int64','read_id':'str'})
	mappeds = mappeds.astype({'chr_rd':'category','start_rd':'int64','end_rd':'int64','read_id':'str'})
	mappeds['idx_map'] = mappeds.index
	mappeds = mappeds.merge(beds_all[['chr_rd','start_rd','end_rd','read_id','str_rd','bc','umi','nb_splices']], how='left')
	mappeds['frag_id'] = mappeds.bc.astype(str) + '/' + mappeds.umi.astype(str)
	#discard all the bc/umi associated to umappeds
	mappeds = mappeds[~mappeds.frag_id.isin(unmappeds)]
	#calculate frag_id number
	# container = pd.DataFrame([[tab[0], len(tab[1])] for tab in mappeds.groupby('frag_id')])
	container = mappeds.groupby('frag_id').size().reset_index()
	container.columns = ['frag_id','frag_id_nb']
	mappeds = mappeds.merge(container, on='frag_id', how='left')

	#exons
	exons_mappeds = pd.DataFrame(bed_gtfs.exons.str.split(';'))
	exons_mappeds['idx'] = exons_mappeds.index
	exons_mappeds = (exons_mappeds.explode('exons'))
	edix = exons_mappeds.idx.tolist()
	exons_mappeds = exons_mappeds.exons.str.split('\t', expand=True)
	exons_mappeds.columns = ['Chromosome','Start','End','gene_name','transcript_name','exon_number']
	exons_mappeds['idx_map'] = edix
	exons_mappeds = exons_mappeds.astype({'Chromosome':'str','Start':'int64','End':'int64','gene_name':'category','transcript_name':'category','exon_number':'int8'})
	exons_mappeds = exons_mappeds.merge(exons, how='left', on = ['Chromosome','Start','End','gene_name','transcript_name','exon_number'])

	#join tables
	mappeds = mappeds.merge(exons_mappeds, on='idx_map', how='left')

	#checkings
	mappeds = mappeds[mappeds.str_rd==mappeds.Strand][['chr_rd','start_rd','end_rd','str_rd','Start','End','StartR','EndR','frag_id','read_id','bc','umi','gene_id','gene_name','transcript_id','transcript_name','exon_id','exon_number','nb_splices','salmon_rlp','frag_id_nb']]
	#discard fragment overlapping out of authorized area
	mappeds['frag_id_transcript'] = mappeds.frag_id.astype('str') + "-" + mappeds.transcript_name.astype('str')
	#coding areas
	to_discard = mappeds[(mappeds.exon_number!=1) & ((mappeds.start_rd < mappeds.Start) | (mappeds.end_rd > mappeds.End))][['frag_id','transcript_name']]
	to_discard = to_discard.frag_id + '-' + to_discard.transcript_name
	mappeds = mappeds[~mappeds.frag_id_transcript.isin(to_discard)]
	#utr areas
	to_discard = mappeds[(mappeds.exon_number==1) & (((mappeds.str_rd=='+') & (mappeds.start_rd < mappeds.Start)) | ((mappeds.str_rd=='-') & (mappeds.end_rd > mappeds.End)))][['frag_id','transcript_name']]
	to_discard = to_discard.frag_id + '-' + to_discard.transcript_name
	mappeds = mappeds[~mappeds.frag_id_transcript.isin(to_discard)]
	#utrs ends
	to_discard = mappeds[(mappeds.exon_number==1) & (((mappeds.str_rd=='+') & (mappeds.end_rd > mappeds.End + MARGIN_END)) | ((mappeds.str_rd=='-') & (mappeds.start_rd < mappeds.Start - MARGIN_END)))][['frag_id','transcript_name']]
	to_discard = to_discard.frag_id + '-' + to_discard.transcript_name
	mappeds = mappeds[~mappeds.frag_id_transcript.isin(to_discard)]


	#frag_id_nb_filtering
	print('checking...')
	mappeds = mappeds.astype({'frag_id':'str','transcript_name':'str','frag_id_nb':'int8'})
	container = mappeds.groupby(['frag_id','transcript_name','frag_id_nb']).size().reset_index()
	container.columns = ['frag_id','transcript_name','frag_id_nb','frag_id_nb_obs']
	container = container[container.frag_id_nb==container.frag_id_nb_obs]
	container['infilter'] = container.frag_id + '-' + container.transcript_name + '-' + container.frag_id_nb.astype(str)
	mappeds['infilter'] = mappeds.frag_id + '-' + mappeds.transcript_name + '-' + mappeds.frag_id_nb.astype(str)
	mappeds = mappeds[mappeds.infilter.isin(container.infilter)]
	del mappeds['infilter']
	# container = [tab[1] for tab in mappeds.groupby(['frag_id','transcript_name']) if (len(tab[1]) == tab[1].frag_id_nb.tolist()[0])]
	# display(mappeds)
	# exit('test')
	# mappeds = pd.concat(container)
	del mappeds['frag_id_transcript']

	#Writing
	mappeds = mappeds.astype({'gene_name':'category', 'frag_id_nb':'int8'})
	mappeds.to_parquet(args.output_path, engine='pyarrow')
