
#import packages
import argparse
import warnings
warnings.filterwarnings("ignore")
import pandas as pd
import gc
import numpy as np
from IPython.display import display
import sys


#Argunent Parser
#***************
parser = argparse.ArgumentParser(description='Map spliced reads on genomic space')
parser.add_argument('IN_path', metavar='INfip', type=str, help='path of overlapped ebed file')
parser.add_argument('transcriptomic_distance_threshold', metavar='trs_dist_value', type=int, default=sys.stdin, help='transcriptomic distance threshold value')
parser.add_argument('output_path', metavar='OUTfop', type=str, help='path of all reads output file')
args = parser.parse_args()

#Functions
#*********
def check_bordering(tab):
	if (np.where((np.diff(x.exon_number)!=1))[0]+1).size == 0:
		return(x)

def checkConsecutive(l):
	n = len(l) - 1
	return (sum(np.diff(sorted(l)) == 1) >= n)


#function to check splicing specificity..
def checking(tab):
	#check2: read_number ~ nsplices
	ltab = len(tab)
	if ltab > 1:
		if ltab == tab.nb_splices.unique()[0]:
			if checkConsecutive(tab.exon_number):
				return(True)
			else:
				return(False)
		else:
			return(False)
	else:
		return(False)


#Variables PARAMS
#****************
TRANSCRIPTOMIC_DISTANCE = args.transcriptomic_distance_threshold
# MARGIN = 1

#files opening
reads = pd.read_parquet(args.IN_path).reset_index(drop=True)

#split spliced and unspliced reads
spliceds = reads[reads.nb_splices>1]
unspliceds = reads[reads.nb_splices==1]

#sorting and casting
spliceds = spliceds.sort_values(['frag_id','read_id','transcript_name','exon_number'])
spliceds = spliceds.astype({'frag_id':'object','read_id':'object'})
#group reads by fragments, read_id, transcript_name
grouped = spliceds[['frag_id','read_id','transcript_name','nb_splices','exon_number']].groupby(['frag_id','read_id','transcript_name'])
spliceds['infilter'] = spliceds.frag_id + '-' + spliceds.read_id + '-' + spliceds.transcript_name
#spliceds = spliceds.astype({'frag_id':'category','read_id':'category'})
save1 = spliceds

#check consecutivity of spliced reads on each transcript for each fragments
print('check consecutivity...')
spliceds_check = (pd.DataFrame([tabs[0][0], tabs[0][1], tabs[0][2], 'splicing_validated']  if checking(tabs[1]) else [tabs[0][0], tabs[0][1], tabs[0][2], 'uncoherent_splicing'] for tabs in grouped))

if len(spliceds_check) != 0:
	spliceds_check.columns = ['frag_id','read_id','trs','splicing_check']
	uncoherent_splicing = spliceds_check[spliceds_check.splicing_check=='uncoherent_splicing']
	spliceds_check = spliceds_check[spliceds_check.splicing_check=='splicing_validated']
	to_remove = uncoherent_splicing.frag_id + '-' + uncoherent_splicing.trs

	#filter validated spliceds into spliceds
	to_keep = spliceds_check.frag_id + '-' + spliceds_check.read_id + '-' + spliceds_check.trs
	spliceds = spliceds[spliceds.infilter.isin(to_keep)]

	#check coordinates bordering of spliced reads
	spliceds = spliceds.astype({'frag_id':'object','transcript_name':'object'})
	uncoherent_bordering = spliceds[(spliceds.start_rd!=spliceds.Start) & (spliceds.end_rd!=spliceds.End)]
	to_remove_bordering = uncoherent_bordering.frag_id + '-' + uncoherent_bordering.transcript_name
	spliceds['infilter'] = spliceds.frag_id + '-' + spliceds.transcript_name
	spliceds = spliceds[~spliceds.infilter.isin(to_remove_bordering)]
	# del spliceds['infilter']

	#filter unspurious transcripts mapping into unspliceds
	unspliceds = unspliceds.astype({'frag_id':'object','read_id':'object'})
	unspliceds['infilter'] = unspliceds.frag_id + '-' + unspliceds.transcript_name
	to_keep = spliceds_check.frag_id + '-' + spliceds_check.trs

	#get from unspliced reads all the fragments not associated with spliced reads
	unspliceds_independant = unspliceds[~unspliceds.frag_id.isin(save1.frag_id)]
	#get from unspliceds reads all the fragments associated with spliced reads validated
	unspliceds_connected = unspliceds[unspliceds.infilter.isin(spliceds.infilter)]

	#merge spliceds and unspliceds
	reads = pd.concat([unspliceds_independant, unspliceds_connected, spliceds])
	del reads['infilter']

else:
	#get from unspliced reads all the fragments not associated with spliced reads
	unspliceds_independant = unspliceds[~unspliceds.frag_id.isin(save1.frag_id)]
	reads = unspliceds_independant




#filter unbalanced frag_id number by transcript
print('unbalance filtering...')
reads = reads.astype({'frag_id':'str', 'transcript_name':'str', 'frag_id_nb':'int8'})
container = reads.groupby(['frag_id','transcript_name','frag_id_nb']).size().reset_index()
container.columns = ['frag_id','transcript_name','frag_id_nb','frag_id_nb_obs']
container = container[container.frag_id_nb==container.frag_id_nb_obs]
container['infilter'] = container.frag_id.astype(str) + '-' + container.transcript_name.astype(str) + '-' + container.frag_id_nb.astype(str)
reads['infilter'] = reads.frag_id.astype(str) + '-' + reads.transcript_name.astype(str) + '-' + reads.frag_id_nb.astype(str)
reads = reads[reads.infilter.isin(container.infilter)]
del reads['infilter']

# reads = pd.concat([tab[1] for tab in reads.groupby(['frag_id','transcript_name']) if len(tab[1]) == tab[1].frag_id_nb.to_list()[0]])

#Calculate reads relative coordinates
#positive strand(+)
tabs = {}
tabs['bed_positive'] = reads[(reads.str_rd == "+")]
tabs['bed_positive']['Start_rdR'] = tabs['bed_positive'].EndR - (tabs['bed_positive'].end_rd - tabs['bed_positive'].Start) + 1
tabs['bed_positive']['End_rdR'] = tabs['bed_positive'].EndR - (tabs['bed_positive'].start_rd - tabs['bed_positive'].Start)
tabs['bed_positive']['dist_END'] = tabs['bed_positive'].Start_rdR
#Negative strand (-)
tabs['bed_negative'] = reads[(reads.str_rd == "-")]
tabs['bed_negative']['Start_rdR'] = tabs['bed_negative'].EndR - (tabs['bed_negative'].End - tabs['bed_negative'].start_rd) + 1
tabs['bed_negative']['End_rdR'] = tabs['bed_negative'].EndR - (tabs['bed_negative'].End - tabs['bed_negative'].end_rd)
tabs['bed_negative']['dist_END'] = tabs['bed_negative'].Start_rdR
# print('concatening..')
reads = pd.concat([tabs['bed_negative'], tabs['bed_positive']])
tabs.clear()

#and now let's filter fragments-transcripts associated to fragment with a distance threshold out of the defined value
to_filter = reads[reads.dist_END > TRANSCRIPTOMIC_DISTANCE]
to_filter = to_filter.frag_id + '-' + to_filter.transcript_name
reads['infilter'] = reads.frag_id + '-' + reads.transcript_name
reads = reads[~reads.infilter.isin(to_filter)]
del reads['infilter']

#delete pcr replicates
reads = reads.drop_duplicates(['chr_rd','start_rd','end_rd','str_rd','Start','End','read_id','transcript_name','frag_id'])
reads = reads.astype({'frag_id':'category','Start_rdR':'int64','End_rdR':'int64','dist_END':'int32','transcript_name':'category'})

#Writing
print("writing...")
reads.to_parquet(args.output_path, engine='pyarrow')
