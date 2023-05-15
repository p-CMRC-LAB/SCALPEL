


import argparse
import pandas as pd
import numpy as np
import vaex as vx
from numba import njit


# Argument Parser
# -----------------
parser = argparse.ArgumentParser(description='Filtering of unmappeds reads and Cleaning')
parser.add_argument('bed_gtf', type=str, help='path of bed_gtf file')
parser.add_argument('exons', type=str, help='path of exons file')
parser.add_argument('tr_distance', type=int, help='transcriptomic distance end threshold value')
parser.add_argument('output_path', type=str, help='path of output bed file')
args = parser.parse_args()

#params
TRANSCRIPTOMIC_DISTANCE = args.tr_distance


# [import file]
# -------------
print("read file...")
mappeds = pd.read_csv(args.bed_gtf, sep="\t", header=None, names=["seqname", "start_rd", "end_rd", "read_id", "strand", "bc", "umi", "nb_splices","exon_id"])
exons = pd.read_csv(args.exons, sep="\t")


# discard unmapped associated fragments (1 ~global mapping)
# -------------------------------------------------------
# add fragID metadata column
mappeds['frag_id'] = mappeds.bc + "-" + mappeds.umi

#delete all the fragment associated to unmapped reads
associated_frags = mappeds[mappeds.exon_id.isnull()].frag_id.tolist()
mappeds = mappeds[~mappeds.frag_id.isin(associated_frags)]
#explode exon_id & drop duplicates
mappeds.exon_id = mappeds.exon_id.str.split(";")
mappeds = mappeds.explode("exon_id").drop_duplicates()


# annotation from GTF via joining by exon_ID
# ----------------------------------------
mappeds = pd.merge(mappeds, exons, on=["seqname","exon_id","strand"], how="left").drop_duplicates().dropna()


# [3] Looks for fragments coherency (1)
# ---------------------------------

#check transcripts out of exons
mappeds = mappeds.astype({"transcript_name":"str"})
mappeds["fidt"] = mappeds.frag_id + "_" + mappeds.transcript_name
associated_frags1 = mappeds[(mappeds.exon_ind!=0) & ((mappeds.start_rd<mappeds.start) | (mappeds.end_rd>mappeds.end))].fidt
associated_frags2 = mappeds[(mappeds.exon_ind==0) & (mappeds.strand=="+") & (mappeds.start_rd<mappeds.start)].fidt
associated_frags3 = mappeds[(mappeds.exon_ind==0) & (mappeds.strand=="-") & (mappeds.end_rd>mappeds.end)].fidt
associated_frags2 = mappeds[(mappeds.exon_ind==0) & (mappeds.strand=="+") & (mappeds.start_rd<mappeds.start)].fidt
associated_frags3 = mappeds[(mappeds.exon_ind==0) & (mappeds.strand=="-") & (mappeds.end_rd>mappeds.end)].fidt
mappeds = mappeds[~mappeds.fidt.isin(associated_frags1)]
mappeds = mappeds[~mappeds.fidt.isin(associated_frags2)]
mappeds = mappeds[~mappeds.fidt.isin(associated_frags3)]

#check fragment coverage
A = mappeds[["transcript_name","frag_id","exon_id"]].drop_duplicates().groupby(["frag_id","transcript_name"]).size().reset_index()
A.columns = ["frag_id","transcript_name","counts"]
A = A.sort_values(["frag_id","counts"])
#get the transcripts for which the frag_id covers the maximum of exons
B = A.drop_duplicates(["frag_id"], keep="last")[["frag_id","counts"]]
B.columns = ["frag_id","max_counts"]
C = A.merge(B, on=["frag_id"], how="left")
C = C[C.counts==C.max_counts]
C["fidt"] = C.frag_id + "_" + C.transcript_name
if len(C) != 0:
	mappeds = mappeds[mappeds.fidt.isin(C.fidt)]



# [4] Discard uncoherent spliced reads
# ------------------------------------
mappeds.read_id = mappeds.read_id.str.replace("/[0-9]","", regex=True)
spliceds = mappeds[mappeds.nb_splices>1]

print("a...")
#discard fragments based on bordering (1)
to_delete = spliceds[(spliceds.start_rd!=spliceds.start) & (spliceds.end_rd!=spliceds.end)].fidt
mappeds = mappeds[~mappeds.fidt.isin(to_delete)]

print("b...")
#discard uncohenrent number spliced reads (1)
A = spliceds.groupby(["fidt","read_id"]).size().reset_index()
A.columns = ["fidt","read_id","counts"]
to_delete = A[A.counts<=1].fidt
mappeds = mappeds[~mappeds.fidt.isin(to_delete)]

print("c...")
#discard uncoherent consecutivity spliced reads (2)
@njit
def group_diff(groups: np.array, values: np.array, lag: int) -> np.array:
	result_exp_mean = np.empty_like(values, dtype=np.float64)
	for i in range(values.shape[0]):
		if groups[i] == groups[i - lag]:
			result_exp_mean[i] = values[i] - values[i - lag]
		else:
			result_exp_mean[i] = np.nan
	return result_exp_mean

A = (spliceds[["fidt","read_id","exon_ind"]]).sort_values(["fidt","read_id","exon_ind"])
groups = A.groupby(["fidt","read_id"]).ngroup().values
values = A.exon_ind.values
A["consecutive"] = group_diff(groups, values, 1)
A["consecutive"] = A.consecutive.fillna(1)
to_delete = A[A.fidt.isin(A[(A.consecutive>1) | (A.consecutive==0)].fidt)].fidt
mappeds = mappeds[~mappeds.fidt.isin(to_delete)]



# [5] Calculate reads relative coordinates
# ------------------------------------
# positive strand(+)
tabs = {}
tabs['bed_positive'] = mappeds[(mappeds.strand == "+")]
tabs['bed_positive']['rel_start_rd'] = tabs['bed_positive'].rel_end + (tabs['bed_positive'].end - tabs['bed_positive'].end_rd)
tabs['bed_positive']['rel_end_rd'] = tabs['bed_positive'].rel_end + (tabs['bed_positive'].end - tabs['bed_positive'].start_rd	)
tabs['bed_positive']['dist_END'] = tabs['bed_positive'].rel_start_rd
#Negative strand (-)
tabs['bed_negative'] = mappeds[(mappeds.strand == "-")]
tabs['bed_negative']['rel_start_rd'] = tabs['bed_negative'].rel_end - (tabs['bed_negative'].end - tabs['bed_negative'].start_rd)
tabs['bed_negative']['rel_end_rd'] = tabs['bed_negative'].rel_end - (tabs['bed_negative'].end - tabs['bed_negative'].end_rd)
tabs['bed_negative']['dist_END'] = tabs['bed_negative'].rel_start_rd
# print('concatening..')
mappeds = pd.concat([tabs['bed_negative'], tabs['bed_positive']])
tabs.clear()


#and now let's filter fragments-transcripts associated to fragment with a distance threshold out of the defined value
mappeds = mappeds[~mappeds.fidt.isin(mappeds[mappeds.dist_END > TRANSCRIPTOMIC_DISTANCE].fidt)].reset_index(drop=True)
mappeds = mappeds.astype({"nb_splices":"str"})
del mappeds["fidt"]
del mappeds["nb_splices"]

#turn Matrix into a vaex table (I/O speed purpose)
mappeds = vx.from_pandas(mappeds)

#Writing
print("writing...")
mappeds.export_hdf5(args.output_path)

# mappeds.to_parquet(args.output_path, engine='fastparquet')
# mappeds.to_csv(args.output_path, sep="\t", header=True, index=False)



