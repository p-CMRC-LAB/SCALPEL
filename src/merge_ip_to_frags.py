
import argparse
import pandas as pd
import numpy as np
import vaex as vx


# Argument Parser
# -----------------
parser = argparse.ArgumentParser(description='Annotate ip positions to fragments')
parser.add_argument('bed_gtf_filtered', metavar='reads', type=str, help='path of bed_gtf filtered file')
parser.add_argument('exons_ips_mappeds', metavar='exons', type=str, help='path of exons associated to mapped_ips file')
parser.add_argument('ip_file', metavar='ipf', type=str, help='path of internal priming ref file')
parser.add_argument('gtf_unique_path', metavar='Efip', type=str, help='path of unique gene file')
parser.add_argument('output_path', metavar='Bmfip', type=str, help='path of output fragment bed file')
parser.add_argument('output_path_ipp', metavar='Bmfip', type=str, help='path of output ip positions bed file')
parser.add_argument('output_bed_unique', metavar='Bmfip', type=str, help='path of output filtered bed unique file')
args = parser.parse_args()

THRESHOLD_SCAN = 200

# [import file] (1)
# -------------
print("read file...")
reads = pd.read_csv(args.bed_gtf_filtered, sep="\t", header=None, names=["chr", "Start_rd", "End_rd", "read_id", "Strand_rd", "bc", "umi", "exon_id","frag_id","Start","End","Strand","gene_name","gene_id","transcript_name","transcript_id","exon_number","StartR","EndR","salmon_rlp","Start_rdR","End_rdR","dist_END"])
exips = pd.read_csv(args.exons_ips_mappeds, sep="\t", header=None, names = ["chr","Start","End","exon_id","Strand","exon_number","transcript_name","ips_index"])
ipfile = pd.read_csv(args.ip_file, sep="\t", header=None, names = ["chr_ip","Start_ip","End_ip","Strand"])

#exons_unique
gtf_unique = pd.read_csv(args.gtf_unique_path, sep="\t", header=None, names=['chr','gene_name','transcript_name'])

# [get ip coords associated to indexs] (2)
# ------------------------------------
exips.ips_index = exips.ips_index.str.split(";")
exips = exips.explode("ips_index")
ipfile["ips_index"] = ipfile.index + 1
ipfile = ipfile.astype({"ips_index":"str"})
exips = exips.merge(ipfile, on=["ips_index","Strand"], how="left")
del exips['ips_index']
del exips['chr_ip']

# Let's discard the IP position in the are of a transcript END [threshold] (3)
# ------------------------------------------------------------
exips["dEND"] = np.where(exips.Strand == "-", exips.Start_ip - exips.Start, exips.End - exips.End_ip)
exips = exips[~((exips.exon_number==1) & (exips.dEND < THRESHOLD_SCAN))]
exips = exips[~exips.Start_ip.isnull()]
exips = exips.sort_values(["transcript_name","exon_number"])
exips = pd.DataFrame([group[1].iloc[0,:].tolist() for group in exips.groupby("transcript_name")])
exips.columns = ["chr","Start","End","exon_id","Strand","exon_number","transcript_name","Start_ip","End_ip","dEND"]


# [Now deal with internal priming filtering] (3)
# -------------------------------------------
reads = reads.merge(exips[["chr","Strand","transcript_name","Start_ip","End_ip","dEND"]], on=["chr","Strand","transcript_name"], how="left")
reads["fidt"] = reads.frag_id + "_" + reads.transcript_name

#reads_noip = reads[reads.Start_ip.isnull()]
reads_ip = reads[~reads.Start_ip.isnull()]
reads_ip_todsc_pos = reads_ip[((reads_ip.Strand=="+") & (reads_ip.Start_rd<=reads_ip.End_ip))]
reads_ip_todsc_neg = reads_ip[((reads_ip.Strand=="-") & (reads_ip.Start_rd>reads_ip.Start_ip))]

if len(reads_ip_todsc_pos) != 0:
	reads = reads[~reads.fidt.isin(reads_ip_todsc_pos.fidt)]

if len(reads_ip_todsc_neg) != 0:
	reads = reads[~reads.fidt.isin(reads_ip_todsc_neg.fidt)]

reads = reads[["chr", "Start_rd", "End_rd", "read_id", "Strand_rd", "bc", "umi", "exon_id","frag_id","Start","End","Strand","gene_name","gene_id","transcript_name","transcript_id","exon_number","StartR","EndR","salmon_rlp","Start_rdR","End_rdR","dist_END"]]
#reads.columns = [['chr_rd','start_rd','end_rd','read_id','str_rd','bc','umi','exon_id','frag_id','Start','End','Strand','gene_name','transcript_name','exon_number','StartR','EndR','salmon_rlp','Start_rdR','End_rdR','dist_END']]
reads = reads.drop_duplicates()


# Filter reads associated to genes with unique transcripts
# -------------------------------------------------------
reads_unique = reads[["gene_name","transcript_name"]].drop_duplicates().groupby("gene_name").count().reset_index(drop=False)
reads_unique = reads_unique[reads_unique.transcript_name == 1]
reads_uniq = reads[reads.gene_name.isin(reads_unique.gene_name)]


# Writing
# -------
vx.from_pandas(reads).export_hdf5(args.output_path)
exips[["Start_ip","End_ip"]].drop_duplicates().to_csv(args.output_path_ipp,sep="\t",index=False)
#unique reads
reads_uniq.to_csv(args.output_bed_unique, sep="\t", doublequote=False, index=None)
