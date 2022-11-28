

#import packages
import argparse
import warnings
warnings.filterwarnings("ignore")
import pandas as pd
import gc
import numpy as np
from IPython.display import display
import sys

#Argument Parser
#**************
parser = argparse.ArgumentParser(description='Filtering of internal priming...')
parser.add_argument('bed_path', metavar='Bfip1', type=str, help='path of fragment file')
parser.add_argument('exons_path', metavar='Efip', type=str, help='path of exon file')
parser.add_argument('ipn_path', metavar='Efip', type=str, help='path of internal priming ref file')
parser.add_argument('gtf_unique_path', metavar='Efip', type=str, help='path of unique gene file')
parser.add_argument('output_bed', metavar='Bmfip', type=str, help='path of output filtered bed file')
parser.add_argument('output_bed_unique', metavar='Bmfip', type=str, help='path of output filtered bed unique file')
parser.add_argument('readid', metavar='Bmfip', type=str, help='path of output readid file')
args = parser.parse_args()


#Defaults variables
IP_THRESHOLD_DISTANCE = 300
ipn_dtypes = {'start_ipn':'int64','end_ipn':'int64','str_rd':'category'}

#Functions
#**********
def fragment_clustering(reads):
	if len(reads)==1:
		return([reads[0][0],reads[0][1],reads[0][2],reads[0][3],reads[0][4],reads[0][5],reads[0][6],[reads[0][7]]])
	else:
		indexes = [x[-1] for x in reads]
		return([reads[0][0],reads[-1][1],reads[0][2],reads[0][3],reads[0][4],reads[0][5],reads[0][6],indexes])

def internal_priming(frag_tab, ip_tab, strand, scope_ip_threshold=300):
	#indexify frag_tab
	frag_tab['frag_index'] = frag_tab.index
	ip_tab['ip_index'] = ip_tab.index
	#initialization
	matched_ips = defaultdict(list)
	j = 0
	#loop on tables
	if strand == '+':
		#reorder tables
		frag_tab = frag_tab.sort_values(by=['End_rd','Start_rd'], ascending=[False,False]).reset_index(drop=True)
		ip_tab = ip_tab.sort_values(by=['End_ipn','Start_ipn'], ascending=[False,False]).reset_index(drop=True)
		#listify tables
		frags = frag_tab[['Start_rd','End_rd','Genomic_End','frag_index']].values.tolist()
		ips = ip_tab[['Start_ipn','End_ipn','ip_index']].values.tolist()
		#loop on the tables
		for ip_index in range(len(ip_tab)):
			current_ip = ips[ip_index]
			#loop on ip table as long than they are superior to the current fragment
			if current_ip[0] - scope_ip_threshold > frags[j][1]:
				continue
			else:
				for frag_index in range(j, len(frag_tab)):
					current_frag = frags[frag_index]
					#loop on frag table as long than they are  superior to the ip_scope (decreasing...)
					if current_frag[1] > current_ip[0]: #+ cross_ip_threshold:
						j += 1
						continue
					else:
						#break loop if frag coord end is inferior to the ip scope
						if current_frag[0] < current_ip[0] - scope_ip_threshold:
							break
						else:
							#check if filling the matching criteria (distance to ip_end / distance to 3'END)
							# if abs(current_ip[1] - current_frag[1]) < abs(current_frag[2] - current_frag[1]) and current_frag[0] < current_ip[1]:
							#New criteria based on bibliographic searches (Article: https://academic.oup.com/nargab/article/4/2/lqac035/6592171) ..consider an internal primig position only if no 3'end in the same scope_ip_threshold...
							if current_ip[1] < current_frag[2] - scope_ip_threshold:
								matched_ips[current_ip[2]].append(current_frag[3])
							else:
								#keep looping on the fragments
								j += 1
								continue
			#stop if no more fragments to scan
			if j >= len(frag_tab):
				break
		#return
		ip_indexs_res = (pd.DataFrame({'ip_index':matched_ips.keys(), 'frag_index': [x[1] for x in matched_ips.items()]})).explode('frag_index')
		frag_tab = frag_tab.merge(ip_indexs_res, how='left')
		frag_tab['ip_index'] = frag_tab['ip_index'].fillna(-1)
		return(frag_tab.sort_values(by=['Start_rd','End_rd','Strand_rd'], ascending=[True,False,True]))
	else:
		#listify tables
		frags = frag_tab[['Start_rd','End_rd','Genomic_End','frag_index']].values.tolist()
		ips = ip_tab[['Start_ipn','End_ipn','ip_index']].values.tolist()
		#loop on the tables
		for ip_index in range(len(ip_tab)):
			current_ip = ips[ip_index]
			#loop on ip table as long than they are inferior to the current fragment
			if current_ip[1] + scope_ip_threshold < frags[j][0]:
				continue
			else:
				for frag_index in range(j, len(frag_tab)):
					current_frag = frags[frag_index]
					#loop on frag table as long than they are  inferior to the ip_scope (decreasing...)
					if current_frag[0] < current_ip[1]: #- cross_ip_threshold:
						j += 1
						continue
					else:
						#break loop if frag coord end is superior to the ip scope
						if current_frag[0] > current_ip[1] + scope_ip_threshold:
							break
						else:
							#check if filling the matching criteria (distance to ip_end / distance to 3'END)
							# if abs(current_frag[0] - current_ip[0]) < abs(current_frag[0] - current_frag[2]) and current_frag[1] > current_ip[0]:
							#New criteria based on bibliographic searches (Article: https://academic.oup.com/nargab/article/4/2/lqac035/6592171) ..consider an internal primig position only if no 3'end in the same scope_ip_threshold...
							if current_ip[0] > current_frag[2] + scope_ip_threshold:
								matched_ips[current_ip[2]].append(current_frag[3])
							else:
								#keep looping on the fragments
								j += 1
								continue
			#stop if no more fragments to scan
			if j >= len(frag_tab):
				break
		#return
		ip_indexs_res = (pd.DataFrame({'ip_index':matched_ips.keys(), 'frag_index': [x[1] for x in matched_ips.items()]})).explode('frag_index')
		frag_tab = frag_tab.merge(ip_indexs_res, how='left')
		frag_tab['ip_index'] = frag_tab['ip_index'].fillna(-1)
		return(frag_tab.sort_values(by=['Start_rd','End_rd','Strand_rd'], ascending=[True,False,True]))


#files Opening
#reads
reads = pd.read_parquet(args.bed_path).reset_index(drop=True)


#ips
ipn = pd.read_csv(args.ipn_path, sep='\*\*\*\*', names = ['ipn','exons'])
ipn = ipn[~ipn.exons.isnull()]
ipn = pd.concat([ipn.ipn.str.split('\t', expand=True), ipn.exons.str.split(';')], axis=1)
ipn = ipn.explode('exons')
# ipn = ipn[[0,1,2,5,'exons']]
ipn.columns = ['chr_ipn','start_ipn','end_ipn','str_ipn','exons']
ipn = pd.concat([ipn[['chr_ipn','start_ipn','end_ipn','str_ipn']], ipn.exons.str.split('\t', expand=True)], axis=1)
ipn.columns = ['chr_ipn','start_ipn','end_ipn','str_rd','chr_ex','Start','End','gene_name','transcript_name','exon_number']
del ipn['chr_ex']
del ipn['exon_number']
del ipn['Start']
del ipn['End']
del ipn['chr_ipn']
ipn = ipn.astype(ipn_dtypes)
ipn = ipn.drop_duplicates()

#exons
exons = pd.read_parquet(args.exons_path)
#exons_unique
gtf_unique = pd.read_csv(args.gtf_unique_path, sep="\t", header=None, names=['gene_name','transcript_name'])


#merge information about last genomic coordinates for each transcript (EXONS / FRAGMENTS)
A = exons[(exons.Strand == "-") & (exons.exon_number == 1)][['transcript_name','Start']].drop_duplicates()
B = exons[(exons.Strand == "+") & (exons.exon_number == 1)][['transcript_name','End']].drop_duplicates()
A.columns = ['transcript_name','Genomic_End']
B.columns = ['transcript_name','Genomic_End']
C = pd.concat([A,B])
reads = reads.merge(C, on = ['transcript_name'], how='left')
del A
del B
del C
del exons
gc.collect()


#Internal priming filtering on each strand of the fragments
reads = reads.astype({'gene_name':'str'})
matched = reads.merge(ipn, on=['str_rd','gene_name','transcript_name'], how='left')

#filtering on matcheds
IP_THR_BACK = 50
matched = matched[((matched.str_rd == '+') & ((matched.start_rd < matched.start_ipn) & (matched.start_rd > matched.start_ipn - IP_THRESHOLD_DISTANCE) & ((matched.Genomic_End > matched.end_ipn + IP_THR_BACK) & (matched.Genomic_End < matched.start_ipn - IP_THRESHOLD_DISTANCE)))) | ((matched.str_rd == '-') & ((matched.end_rd > matched.end_ipn) & (matched.end_rd < matched.end_ipn + IP_THRESHOLD_DISTANCE) & ((matched.Genomic_End > matched.end_ipn + IP_THRESHOLD_DISTANCE) & (matched.Genomic_End < matched.start_ipn - IP_THR_BACK))))]

reads['infilter'] = reads.frag_id.astype(str) + '-' + reads.transcript_name.astype(str)
reads = reads.astype({'transcript_name':'category'})
if len(matched) != 0:
	matched = matched.frag_id.astype(str) + '-' + matched.transcript_name.astype(str)
	reads = reads[~reads.infilter.isin(matched)]
else:
	pass

#filter out reads mapping out of coordinates in coding area
out_reads = reads[(reads.exon_number!=1) & ((reads.start_rd<reads.Start-1) | (reads.end_rd>reads.End+1))]
reads = reads[~reads.infilter.isin(out_reads.infilter)]
#filter out reads mapping out of coordinates in utr area
out_reads = reads[(reads.exon_number==1) & (reads.str_rd=='+') & ((reads.start_rd<reads.Start-1))]
reads = reads[~reads.infilter.isin(out_reads.infilter)]
out_reads = reads[(reads.exon_number==1) & (reads.str_rd=='-') & ((reads.end_rd>reads.End+1))]
reads = reads[~reads.infilter.isin(out_reads.infilter)]
del reads['infilter']
del reads['Genomic_End']

#generate reads mapped on unique fragments
reads_uniq = reads[reads.gene_name.isin(gtf_unique.gene_name.tolist())]


#writing
#read_id/bc/umi
(reads[['read_id','bc','umi']].drop_duplicates()).to_parquet(args.readid, engine='pyarrow')
reads = reads[['bc','gene_name','gene_id','transcript_name','transcript_id','salmon_rlp','umi','frag_id','dist_END']]
#unique reads
reads_uniq.to_csv(args.output_bed_unique, sep="\t", doublequote=False, index=None)
#reads
reads.to_parquet(args.output_bed, engine='pyarrow')




#
             
