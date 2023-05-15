

import warnings
warnings.filterwarnings("ignore")
import argparse
import pandas as pd
import numpy as np
from functools import reduce
import operator


#Argunent Parser
#**************
parser = argparse.ArgumentParser(description='Processing of GTF file')
parser.add_argument('gtf_file', type=str, help='path of gtf file splitted by chromosome')
parser.add_argument('transcriptomic_distance_threshold', type=int, help='transcriptomic distance threshold value')
parser.add_argument('transcript_end_distance_threshold', type=int, help='transcriptomic distance end threshold value')
parser.add_argument('exon_file', type=str, help='path of the exon entries output bed file')
parser.add_argument('exon_unique_file', type=str, help='path of the unique gene exon entries output bed file')
args = parser.parse_args()

#params
TRANSCRIPTOMIC_DISTANCE = int(args.transcriptomic_distance_threshold)
THRESHOLD_TRANSCRIPT_END = int(args.transcript_end_distance_threshold)

#custom functions
def convert_to_relative_coords(df):
    A = (df.end - df.start).tolist()
    OUT_start = []
    OUT_end = []
    strand = df.strand.tolist()[0]
    if len(df) > 1:
        for ind, intv in enumerate(A):
            if ind == 0:
                if strand == "-":
                    OUT_start.append(ind)
                    OUT_end.append(ind+intv)
                else:
                    OUT_start.append(ind+intv)
                    OUT_end.append(ind)
            else:
                RES = [A[ind-1]+1, A[ind-1]+intv+1]
                if strand == "-":
                    OUT_start.append(RES[0])
                    OUT_end.append(RES[1])
                else:
                    OUT_start.append(RES[1])
                    OUT_end.append(RES[0])
                #update
                A[ind] = RES[1]
        df["rel_start"] = OUT_start
        df["rel_end"] = OUT_end
    else:
        if strand == "-":
            df["rel_start"] = 0
            df["rel_end"] = df.end - df.start
        else:
            df["rel_start"] = df.end - df.start
            df["rel_end"] = 0

    #add annotation about isoform genomic end coordinate
    if df.strand.tolist()[0] == "-":
        df["end_isoform"] = df.start.tolist()[0]
        df["size_isoform"] = df.rel_end.max()
    else:
        df["end_isoform"] = df.end.tolist()[0]
        df["size_isoform"] = df.rel_start.max()

    #add an index to to the isoform
    df['exon_ind'] = list(range(len(df)))

    #return
    return(df)


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


def get_similar_transcripts(current_g, end_thr, dist_thr):
	"""
	"""
	strand = current_g.strand.tolist()[0]
	transcript_similars_todelete = []
	similar_transcripts = []
	END_THRESHOLD = end_thr
	TRANSCRIPTOMIC_DISTANCE_THRESHOLD = dist_thr

	#loop on exons
	for ind in range(len(current_g.exon_ind.unique())):
		current_ex = current_g[current_g.exon_ind == ind]

		#case0: last exon
		if ind == 0:
			#get similars ends
			if strand == "-":
				cluster_starts = dict(enumerate(grouper(np.sort(current_ex.start), END_THRESHOLD),1))
				similar_trends = [current_ex[current_ex.start.isin(coords)].transcript_name.tolist() for coords in cluster_starts.values()]
				for trs in similar_trends:
					if len(trs) > 1:
						current_ex_sims = current_ex[current_ex.transcript_name.isin(trs)]
						current_ex_sims.end = np.where(((current_ex_sims['rel_end'] > TRANSCRIPTOMIC_DISTANCE_THRESHOLD)), 0, current_ex_sims['end'])
						for hit in current_ex_sims.groupby("end"):
							current_hit = hit[1]
							#check table length and threshold distance
							if (len(current_hit) > 1) and current_hit.rel_end.tolist()[0] > TRANSCRIPTOMIC_DISTANCE_THRESHOLD:
								transcript_similars_todelete.append(current_hit[current_hit.size_isoform != current_hit.size_isoform.max()].transcript_name.tolist())
							elif len(current_hit) > 1 and current_hit.rel_end.tolist()[0] <= TRANSCRIPTOMIC_DISTANCE_THRESHOLD:
								similar_transcripts.append(current_hit.transcript_name.tolist())
							else:
								pass
			else:
				cluster_ends = dict(enumerate(grouper(np.sort(current_ex.end), END_THRESHOLD),1))
				similar_trends = [current_ex[current_ex.end.isin(coords)].transcript_name.tolist() for coords in cluster_ends.values()]
				for trs in similar_trends:
					if len(trs) > 1:
						current_ex_sims = current_ex[current_ex.transcript_name.isin(trs)]
						current_ex_sims.start = np.where((current_ex_sims['rel_start'] > TRANSCRIPTOMIC_DISTANCE_THRESHOLD), 0, current_ex_sims['start'])
						for hit in current_ex_sims.groupby("start"):
							current_hit = hit[1]
							#check table length and threshold distance
							if (len(current_hit) > 1) and current_hit.rel_start.tolist()[0] > TRANSCRIPTOMIC_DISTANCE_THRESHOLD:
								transcript_similars_todelete.append(current_hit[current_hit.size_isoform != current_hit.size_isoform.max()].transcript_name.tolist())
							elif len(current_hit) > 1 and current_hit.rel_start.tolist()[0] <= TRANSCRIPTOMIC_DISTANCE_THRESHOLD:
								similar_transcripts.append(current_hit.transcript_name.tolist())
							else:
								pass
		else:
			#let's work on similars transcript if there is
			updated_similar_transcripts = []
			for ind1, trs in enumerate(similar_transcripts):
				if len(trs) > 1:
					current_hit = current_ex[(current_ex.exon_ind == ind) & (current_ex.transcript_name.isin(trs))]
					#check table and threshold distance
					if strand == "-":
						current_hit.end = np.where((current_hit['rel_end'] > TRANSCRIPTOMIC_DISTANCE_THRESHOLD), current_hit['start'] + TRANSCRIPTOMIC_DISTANCE_THRESHOLD, current_hit['end'])
						for hit1 in current_hit.groupby(["start","end"]):
							current_hit1 = hit1[1]
							if len(current_hit1) > 1:
								#check table length and threshold distance
								if current_hit1.rel_end.tolist()[0] > TRANSCRIPTOMIC_DISTANCE_THRESHOLD:
									transcript_similars_todelete.append(current_hit1[current_hit1.size_isoform != current_hit1.size_isoform.max()].transcript_name.tolist())
								else:
									updated_similar_transcripts.append(current_hit1.transcript_name.tolist())
							else:
								pass
					else:
						current_hit.start = np.where((current_hit['rel_start'] > TRANSCRIPTOMIC_DISTANCE_THRESHOLD), current_hit['end'] - TRANSCRIPTOMIC_DISTANCE_THRESHOLD, current_hit['start'])
						for hit1 in current_hit.groupby(["start","end"]):
							current_hit1 = hit1[1]
							if len(current_hit1) > 1:
								#check table length and threshold distance
								if current_hit1.rel_end.tolist()[0] > TRANSCRIPTOMIC_DISTANCE_THRESHOLD:
									transcript_similars_todelete.append(current_hit1[current_hit1.size_isoform != current_hit1.size_isoform.max()].transcript_name.tolist())
								else:
									updated_similar_transcripts.append(current_hit1.transcript_name.tolist())
							else:
								pass
				else:
					pass
			#remove None occurences
			similar_transcripts = updated_similar_transcripts
			#break the main loop if no more similar transcript
			if len(similar_transcripts) == 0:
				if len(transcript_similars_todelete) == 0:
					return None
				else:
					return(reduce(operator.concat, transcript_similars_todelete))



# ****
# Main
# ****

#File Opening
print("file opening...")
gtf = pd.read_csv(args.gtf_file, sep="\t")

#filter exon & UTR features
gtf_exs = gtf[gtf.feature=="exon"]
gtf_utrs = gtf[gtf.feature=="UTR"]

#add UTR annotation in exon gtf
print("add UTR annotation in gtf file...")
gtf = pd.merge(gtf_exs,gtf_utrs[['exon_id','feature']], on="exon_id", how='left')
gtf.feature_y = gtf.feature_y.fillna("")
gtf['features'] = gtf.feature_x + ";" + gtf.feature_y
gtf = gtf.dropna().drop_duplicates().drop(["feature_x", "feature_y"], axis=1)

#Cast integer columns
gtf = gtf.astype({"start":"int32","end":"int32","exon_number":"int8"})

#add relative coords
print("add relative coords...")
gtf = gtf.sort_values(["transcript_name","exon_number"], ascending=[True,False]).reset_index(drop=True)
gtf = gtf.groupby("transcript_name").apply(lambda x: convert_to_relative_coords(x))

#get similars transcripts
print("get similar transcripts")
sim_transcripts = [get_similar_transcripts(gtab[1], THRESHOLD_TRANSCRIPT_END, TRANSCRIPTOMIC_DISTANCE) for gtab in gtf.groupby(["gene_name"])]
sim_transcripts = [trs for trs in sim_transcripts if trs != None]
if len(sim_transcripts) != 0:
	sim_transcripts = reduce(operator.concat, sim_transcripts)
	gtf = gtf[~gtf.transcript_name.isin(sim_transcripts)]

#Get unique genes associated to chr
gtf_unique = (gtf.groupby(["gene_name","transcript_name"]).size()).groupby("gene_name").size().reset_index()
gtf_unique.columns = ['gene_name','counts']
gtf_unique = gtf_unique[gtf_unique.counts==1]

#writing
print("writing...")
#exons
gtf.to_csv(args.exon_file,sep="\t",index=False)
#unique gene
gtf_unique.to_csv(args.exon_unique_file, sep="\t", index=False)

