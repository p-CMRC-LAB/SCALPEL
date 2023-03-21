

import argparse
import pandas as pd
import numpy as np
from numba import njit
import sys


# Argument Parser
# -----------------
parser = argparse.ArgumentParser(description='Filtering of unmappeds reads and Cleaning')
parser.add_argument('bed_gtf', metavar='reads', type=str, help='path of bed_gtf file')
parser.add_argument('exons', metavar='exons', type=str, help='path of exons file')
parser.add_argument('cross_threshold', type=int, default=sys.stdin, help='transcriptomic distance end threshold value')
parser.add_argument('output_path', metavar='Bmfip', type=str, help='path of output bed file')
args = parser.parse_args()

CROSS_THRESHOLD = 50

# [import file]
# -------------
print("read file...")
bed_gtf = {
	"file": pd.read_csv(args.bed_gtf, sep="\t", names=["chr", "Start_rd", "End_rd", "read_id", "Strand_rd", "bc", "umi", "nb_splices","exon_id"], dtype={"exon_id": "category"})
}
exons = pd.read_parquet(args.exons)

# add fragID metadata column
bed_gtf["file"]['frag_id'] = bed_gtf["file"].bc + "-" + bed_gtf["file"].umi

print("1")
# [1] Extract mapped and unmapped reads (1) and discard unmapped associated fragments (1 ~global mapping)
# ----------------------------------------------------------------------------------------------------------
mappeds = bed_gtf["file"][~bed_gtf["file"].exon_id.isnull()]
unmappeds = bed_gtf["file"][bed_gtf["file"].exon_id.isnull()][['frag_id']]
bed_gtf.clear()
# discarding
mappeds = mappeds[~mappeds.frag_id.isin(unmappeds.frag_id)]


print("2")
# [2] split the table by exon_id and add annotation from GTF via joining by exon_ID
# ----------------------------------------------------------------------------------
# splitting
mappeds.exon_id = mappeds.exon_id.str.split(";")
mappeds = mappeds.explode("exon_id").drop_duplicates()
# add transcript_name metadata
mappeds = mappeds.merge(exons[["exon_id","Start","End","Strand","gene_name","gene_id","transcript_name","transcript_id","exon_number","StartR","EndR","salmon_rlp"]], on = "exon_id", how="left")
mappeds = mappeds.astype({"transcript_name":"str"})
mappeds = mappeds[mappeds.Strand_rd==mappeds.Strand]
#remove reads which overlap out of exons
mappeds["fidt"] = mappeds.frag_id + "_" + mappeds.transcript_name
target = mappeds[(mappeds.exon_number!=1) & ((mappeds.Start_rd<mappeds.Start) | (mappeds.End_rd > mappeds.End))]
mappeds = mappeds[~mappeds.fidt.isin(target.fidt)]
mappeds_pos_target = mappeds[(mappeds.exon_number==1) & (mappeds.Strand_rd == "+") & ((mappeds.Start_rd < mappeds.Start) | (mappeds.End_rd > (mappeds.End + CROSS_THRESHOLD)))]
mappeds_neg_target = mappeds[(mappeds.exon_number==1) & (mappeds.Strand_rd == "-") & ((mappeds.End_rd > mappeds.End) | (mappeds.Start_rd < (mappeds.Start - CROSS_THRESHOLD)))]
mappeds = mappeds[~mappeds.fidt.isin(mappeds_pos_target.fidt)]
mappeds = mappeds[~mappeds.fidt.isin(mappeds_neg_target.fidt)]


print("3")
# [3] Looks for fragments coherency (1)
# ---------------------------------

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


print("4")
# [4] Discard uncoherent spliced reads
# ------------------------------------
mappeds.read_id = mappeds.read_id.str.replace("/[0-9]","", regex=True)
spliceds = mappeds[mappeds.nb_splices>1]

print("a...")
#discard fragments based on bordering (1)
to_delete = spliceds[(spliceds.Start_rd!=spliceds.Start) & (spliceds.End_rd!=spliceds.End)].fidt
mappeds = mappeds[~mappeds.fidt.isin(to_delete)]

print("b...")
spliceds = mappeds[mappeds.nb_splices>1]
#discard uncohenrent number spliced reads (1)
A = spliceds.groupby(["fidt","read_id"]).size().reset_index()
A.columns = ["fidt","read_id","counts"]
to_delete = A[A.counts<=1].fidt
mappeds = mappeds[~mappeds.fidt.isin(to_delete)]

print("c...")
spliceds = mappeds[mappeds.nb_splices>1]
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

A = (spliceds[["fidt","read_id","exon_number"]]).sort_values(["fidt","read_id","exon_number"])
groups = A.groupby(["fidt","read_id"]).ngroup().values
values = A.exon_number.values
A["consecutive"] = group_diff(groups, values, 1)
A["consecutive"] = A.consecutive.fillna(1)
to_delete = A[A.fidt.isin(A[(A.consecutive>1) | (A.consecutive==0)].fidt)].fidt
mappeds = mappeds[~mappeds.fidt.isin(to_delete)]



print("5")
# Calculate reads relative coordinates
# ------------------------------------
# positive strand(+)
tabs = {}
tabs['bed_positive'] = mappeds[(mappeds.Strand_rd == "+")]
tabs['bed_positive']['Start_rdR'] = tabs['bed_positive'].EndR - (tabs['bed_positive'].End_rd - tabs['bed_positive'].Start) + 1
tabs['bed_positive']['End_rdR'] = tabs['bed_positive'].EndR - (tabs['bed_positive'].Start_rd - tabs['bed_positive'].Start)
tabs['bed_positive']['dist_END'] = tabs['bed_positive'].Start_rdR
#Negative strand (-)
tabs['bed_negative'] = mappeds[(mappeds.Strand_rd == "-")]
tabs['bed_negative']['Start_rdR'] = tabs['bed_negative'].EndR - (tabs['bed_negative'].End - tabs['bed_negative'].Start_rd) + 1
tabs['bed_negative']['End_rdR'] = tabs['bed_negative'].EndR - (tabs['bed_negative'].End - tabs['bed_negative'].End_rd)
tabs['bed_negative']['dist_END'] = tabs['bed_negative'].Start_rdR
# print('concatening..')
mappeds = pd.concat([tabs['bed_negative'], tabs['bed_positive']])
tabs.clear()

#and now let's filter fragments-transcripts associated to fragment with a distance threshold out of the defined value
TRANSCRIPTOMIC_DISTANCE = 600
mappeds = mappeds[~mappeds.fidt.isin(mappeds[mappeds.dist_END > TRANSCRIPTOMIC_DISTANCE].fidt)]
mappeds = mappeds.astype({"nb_splices":"str"})
del mappeds["fidt"]
del mappeds["nb_splices"]

#Writing
print("writing...")
#mappeds.to_parquet(args.output_path, engine='fastparquet')
mappeds.to_csv(args.output_path, sep="\t", header=False, index=False)
