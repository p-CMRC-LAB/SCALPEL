#/usr/bin/python3


# This script have for function to prefilter the internal priming entries (into the exons coordinates scope) and
# remove those close to transcript END coordinate base on a threshold distance value



#libraries importations
#**********************
import warnings
warnings.filterwarnings("ignore")
import argparse
import pyranges as pr
import pandas as pd
import sys
import numpy as np
import functools
import operator
from collections import defaultdict
from IPython.display import display

#Functions
#*********

def relative_coords(tab):
	#reverse exon number
	tab.loc[:,'exon_number'] = tab.exon_number.tolist()[::-1]
	#initialization
	starts = tab.Start.tolist()
	ends = tab.End.tolist()
	new_starts = [None] * len(starts)
	new_ends = [None] * len(ends)
	strand = tab.Strand.tolist()[0]
	#loop
	if strand == "-":
		for ind in range(len(starts)):
			if ind == 0:
				#initialize first_coords
				new_starts[0] = 0
				new_ends[0] = ends[0] - starts[0] -1
			else:
				new_starts[ind] = new_ends[ind-1] +1
				new_ends[ind] = new_starts[ind] + (ends[ind] - starts[ind] -1)
		rl_length = new_ends[-1]
	else:
		for i, ind in enumerate(range(len(starts))[::-1]):
			if i == 0:
				#initialize first_coords
				new_starts[ind] = 0
				new_ends[ind] = ends[ind] - starts[ind] -1
			else:
				new_starts[ind] = new_ends[ind+1] +1
				new_ends[ind] = new_starts[ind] + (ends[ind] - starts[ind] -1)
		rl_length = new_ends[0]
	#return:
	tab.loc[:,'StartR'] = new_starts
	tab.loc[:,'EndR'] = new_ends
	tab.loc[:,'rl_length'] = rl_length
	return(tab)

def grouper(iterable, threshold):
	prev = None
	group = []
	for item in iterable:
		if prev is None or item - prev <= threshold:
			group.append(item)
		else:
			yield group
			group = [item]
		prev = item
	if group:
		yield group

def find_similar_ends(tab, strand, exon_thr):
	if strand == "-":
		#get cluster of start coordinates based on range
		cluster_starts = dict(enumerate(grouper(np.sort(tab.Start), exon_thr),1))
		return([tab[tab.Start.isin(start)].transcript_name.to_list() for start in cluster_starts.values() if len(start)>1])
	else:
		#get clusters of end coordinate based on range
		cluster_ends = dict(enumerate(grouper(np.sort(tab.End), exon_thr),1))
		return([tab[tab.End.isin(end)].transcript_name.to_list() for end in cluster_ends.values() if len(end)>1])

def similar_last_exon(tab, trs_distance_thr, end_exon_thr):
	#Get gene strand
	strand = tab.Strand.tolist()[0]
	outs = []
	if strand == "-":
		#Finf the sets of transcripts with similar ends
		results = find_similar_ends(tab, strand = "-", exon_thr=end_exon_thr)
		if len(results) != 0:
			#Loop on each group of transcript with similar end
			for res in results:
				current_tab = tab[tab.transcript_name.isin(res)]
				#check coord_limit (And set a coord limit if coordinates is superior to the input threshold)
				coord_limit = current_tab.End.tolist()[0] - (current_tab.EndR.to_list()[0] - trs_distance_thr)
				current_tab['End'] = np.where((current_tab['EndR'] > trs_distance_thr), coord_limit, current_tab['End'])
				#group the transcripts with similar 3' end (within threshold) and exact 5'end coordinates
				groups = current_tab.groupby('End')
				for group in groups:
					curr_group = group[1]
					# for each group of similar transcripts if several transcripts
					if len(curr_group) > 1:
						#keep information about the distance between the transcript 3'end VS the shortest (further tasks...)
						# If we reach at this step the transcriptomic distance threshold defined in input... no need to keep looping on the other exons (CHOOSE ONE TRANSCRIPT!)
						if curr_group.EndR.tolist()[0] > trs_distance_thr:
							#Choose one transcript based on the salmon quantification NumReads
							trs_to_keep = (curr_group[curr_group.NumReads == max(curr_group.NumReads)].transcript_name.to_list())[0]
							outs.append([curr_group[curr_group.transcript_name != trs_to_keep].transcript_name.to_list(), "to_discard"])
						# Else indicate in the output then a looping on the other exons is required
						else:
							outs.append([group[1].transcript_name, "to_loop_on"])
		return(outs)
	else:
		results = find_similar_ends(tab, strand = "+", exon_thr=end_exon_thr)
		if len(results) != 0:
			for res in results:
				current_tab = tab[tab.transcript_name.isin(res)]
				#check coord limit
				coord_limit = current_tab.Start.tolist()[0] + (current_tab.EndR.tolist()[0] + trs_distance_thr)
				current_tab['Start'] = np.where((current_tab['EndR'] > trs_distance_thr), coord_limit, current_tab['Start'])
				groups = current_tab.groupby('Start')
				for group in groups:
					curr_group = group[1]
					if len(curr_group) > 1:
						#keep information about the distance between the transcript 3'end VS the shortest (further tasks...)
						# If we reach at this step the transcriptomic distance threshold defined in input... no need to keep looping on the other exons (CHOOSE ONE TRANSCRIPT!)
						if curr_group.EndR.tolist()[0] > trs_distance_thr:
							#Choose one transcript based on the salmon quantification NumReads
							trs_to_keep = (curr_group[curr_group.NumReads == max(curr_group.NumReads)].transcript_name.to_list())[0]
							outs.append([curr_group[curr_group.transcript_name != trs_to_keep].transcript_name.to_list(), "to_discard"])
						# Else indicate in the output then a looping on the other exons is required
						else:
							outs.append([group[1].transcript_name, "to_loop_on"])
		return(outs)

def similar_exon(tab, trs_distance_thr):
	#get strand
	strand = tab.Strand.tolist()[0]
	outs = defaultdict(list)
	#negative strand case
	if strand == "-":
		#1/group the table by Start coordinate
		group_starts = tab.groupby('Start')
		#2/loop on transcripts with same starts
		for group in group_starts:
			#if several transcripts with the same start coordinate
			if len(group[1]) > 1:
				#check the coord limit end criteria
				current_tab = group[1]
				coord_end_limit = current_tab.End.tolist()[0] - (current_tab.EndR.tolist()[0] - trs_distance_thr)
				current_tab['End'] = np.where((current_tab['EndR'] > trs_distance_thr), coord_end_limit, current_tab['End'])
				#group the transcript with same ends
				group_ends = current_tab.groupby('End')
				for group2 in group_ends:
					current_tab2 = group2[1]
					#check if several transcripts with same start and end coordinates
					if len(current_tab2) == 1:
						continue
					else:
						# If we reach at this step the transcriptomic distance threshold defined in input... no need to keep looping on the other exons (CHOOSE ONE TRANSCRIPT! / DISCARD THE OTHERS)
						if current_tab2.EndR.tolist()[0] > trs_distance_thr:
							#Choose one transcript based on the salmon quantification NumReads parameters if NumReads != None
							trs_to_keep = (current_tab2[current_tab2.NumReads == max(current_tab2.NumReads)].transcript_name.to_list())[0]
							outs['to_discard'].append(current_tab2[current_tab2.transcript_name != trs_to_keep].transcript_name.to_list())
							# outs['to_discard'].append((current_tab2.sort_values(by=['rl_length'], ascending=[False])).transcript_name.tolist()[1:])
						else:
							outs['to_loop_on'].append(current_tab2.transcript_name.tolist())
			else:
				#if transcript with an unique start coordinate
				continue
	else:
		#1/group the table by End coordinate
		group_ends = tab.groupby('End')
		#2/loop on transcripts with similar ends
		for group in group_ends:
			#if several transcripts with the same end coordinate
			if len(group[1]) > 1:
				#check the coord limit start criteria
				current_tab = group[1]
				coord_start_limit = current_tab.Start.tolist()[0] + (current_tab.EndR.tolist()[0] - trs_distance_thr)
				current_tab['Start'] = np.where((current_tab['EndR'] > trs_distance_thr), coord_start_limit, current_tab['Start'])
				#group the transcript with similar start
				group_starts = current_tab.groupby('Start')
				for group2 in group_starts:
					current_tab2 = group2[1]
					#check if several transcripts with same start and end coordinates
					if len(current_tab2) == 1:
						continue
					else:
						# If we reach at this step the transcriptomic distance threshold defined in input... no need to keep looping on the other exons (CHOOSE ONE TRANSCRIPT! / DISCARD THE OTHERS)
						if current_tab2.EndR.tolist()[0] > trs_distance_thr:
							#Choose one transcript based on the salmon quantification NumReads parameters if NumReads != None
							trs_to_keep = (current_tab2[current_tab2.NumReads == max(current_tab2.NumReads)].transcript_name.to_list())[0]
							outs['to_discard'].append(current_tab2[current_tab2.transcript_name != trs_to_keep].transcript_name.to_list())
							# outs['to_discard'].append((current_tab2.sort_values(by=['rl_length'], ascending=[False])).transcript_name.tolist()[1:])
						else:
							outs['to_loop_on'].append(current_tab2.transcript_name.tolist())
			else:
				#if transcript with an unique end coordinate
				continue
	return(outs)

def remove_sim_exons(tab, distance_thr, exon_thr):
	#print gene name
	# display(tab.gene_name.tolist()[0])
	#function to remove similars exons in the threshold scope
	tab.loc[:,'exon_number'] = pd.to_numeric(tab['exon_number'])
	exons_number = list(set(tab.exon_number.tolist()))
	exons_number.sort()
	to_discard = []

	max_exnumber = max(tab.exon_number.tolist())

	#Here, check than a gene does contains a last exon. An error can occur if the gene does not have a last exon (exon number == 1)
	if len(tab[tab.exon_number==1]) == 0:
		warnings.warn("Error: Gene " + tab.gene_name.tolist()[0] + " does not contain a last exon [exon_number==1]. Check input format of GTF file !!!")
		return(None)

	# Get list of similar transcripts in last exon
	results = similar_last_exon(tab[tab.exon_number == 1], trs_distance_thr=distance_thr, end_exon_thr=exon_thr)
	if len(results) == 0:
		return(tab)
	else:
		# Loop on transcripts with similar transcriptomic ends with 2 option: (1- keep looking on the next exon) / (2- Pick one transcripts among the similars)
		for result in results:
			to_keeps = []
			ex_number = 2
			#CASE1: Discard transcripts among those with similar ends in the last exon and threshold distance reached
			if result[1] == 'to_discard':
				tab = tab[~tab.transcript_name.isin(result[0])]
			#CASE2: Loop into the next exons for transcripts with similar ends and threshold distance not reached
			elif result[1] == 'to_loop_on':
				exons_similar_transcripts_tab = tab[tab.transcript_name.isin(result[0])]
				checking = True
				# Looping STOP condition: no more result suggesting to keep on the loop
				while ex_number <= max_exnumber:
					# Get the exon table corresponding to the current exon number to check
					exons_similar_transcripts_exnumber_tab = exons_similar_transcripts_tab[exons_similar_transcripts_tab.exon_number == ex_number]
					#if no exons at all
					if len(exons_similar_transcripts_exnumber_tab) == 0:
						break
					else:
						# Get similar exons with 2 options: (1- keep looking on the next exons) / (2- Discard those with similar ends from the input table)
						results_sims = similar_exon(exons_similar_transcripts_exnumber_tab, trs_distance_thr = distance_thr)
						#if no similar exons at all
						if len(results_sims) == 0:
							break
						else:
							#check those to discard from the input table
							if len(results_sims['to_discard']) != 0:
								to_discs = functools.reduce(operator.iconcat, results_sims['to_discard'], [])
								tab = tab[~tab.transcript_name.isin(to_discs)]
							else:
								pass

							#check those for which we have to loop on
							if len(results_sims['to_loop_on']) != 0:
								to_loops = functools.reduce(operator.iconcat, results_sims['to_loop_on'], [])
								exons_similar_transcripts_tab = exons_similar_transcripts_tab[exons_similar_transcripts_tab.transcript_name.isin(to_loops)]
							else:
								break
						#increment exon number if not maximal exon reached
						ex_number += 1
			else:
				raise NameError("Error Function in processing similar_last_exon for the gene: " + tab.gene_name.tolist()[0])
		return(tab)




#Argunent Parser
#**************
parser = argparse.ArgumentParser(description='Processing of GTF file')
parser.add_argument('gtf_file', metavar='GTF_file_input', type=str, default=sys.stdin, help='path of gtf file splitted by chromosome')
parser.add_argument('transcriptomic_distance_threshold', type=int, default=sys.stdin, help='transcriptomic distance threshold value')
parser.add_argument('transcript_end_distance_threshold', type=int, default=sys.stdin, help='transcriptomic distance end threshold value')
parser.add_argument('exon_file', type=str, default=sys.stdout, help='path of the exon entries output bed file')
parser.add_argument('exon_unique_file', type=str, default=sys.stdout, help='path of the unique gene exon entries output bed file')
parser.add_argument('exon_bedmap', type=str, default=sys.stdout, help='path of the exon coords for bedmap output bed file')
args = parser.parse_args()

#Variables
#*********
TRANSCRIPTOMIC_DISTANCE = int(args.transcriptomic_distance_threshold)
THRESHOLD_TRANSCRIPT_END = int(args.transcript_end_distance_threshold)


#Files Opening
#*************a
gtf = pd.read_csv(args.gtf_file, sep="\t")

#Main
#****

#Calculate relative coordinates
print("Calculate relative coordinates...")
gtf = gtf.groupby('transcript_name').apply(lambda x: relative_coords(x))


#Discard all the transcripts with similar ends for each gene
print("Discard all the transcripts with similar ends for each gene...")
gtf = gtf.groupby('gene_name').apply(lambda x: remove_sim_exons(x, distance_thr=TRANSCRIPTOMIC_DISTANCE, exon_thr=THRESHOLD_TRANSCRIPT_END))
gtf = gtf.dropna()

#Sorting Gtf file
print("Sorting Gtf file...")
gtf = gtf.sort_values(by=['Chromosome','Start','End']).reset_index(drop=True)
ex_coltypes = {'Chromosome':'str','Start':'int64','End':'int64','Strand':'str','StartR':'int32','EndR':'int32','gene_id':'category','gene_name':'category','transcript_id':'category','transcript_name':'category','exon_id':'category','exon_number':'int8','gene_type':'category','transcript_type':'category','salmon_rlp':'float'}
gtf = gtf.astype(ex_coltypes)

#Get unique genes associated to chr
print("Get unique genes associated to chr...")
gtf_unique = gtf[['Chromosome','gene_name','transcript_name']].sort_values(by=['Chromosome','gene_name','transcript_name'])
gtf_unique = gtf_unique[~gtf_unique.duplicated(subset='gene_name', keep=False)]


#Writing
print("Writing...")
gtf.to_parquet(args.exon_file, engine='pyarrow')
gtf_unique.to_csv(args.exon_unique_file, sep='\t', doublequote=False, index=False)
gtf[['Chromosome','Start','End','exon_id']].to_csv(args.exon_bedmap, sep='\t', header=False, index=False)
gtf['Chromosome'] = gtf.Chromosome.str.replace('chr','')
