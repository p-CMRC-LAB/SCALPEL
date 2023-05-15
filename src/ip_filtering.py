


#import libs
import argparse
import pandas as pd
import numpy as np
import vaex as vx


# Argument Parser
# -----------------
parser = argparse.ArgumentParser(description='Filter internal priming position')
parser.add_argument('bed_gtf_filtered', type=str, help='path of bed_gtf_filtered file')
parser.add_argument('ip_file', type=str, help='path of internal priming ref file')
parser.add_argument('output_path', type=str, help='path of output file for ip filtered reads')
parser.add_argument('output_bed_unique', type=str, help='path of output file for read mapped in unique transcript gene')
parser.add_argument('read_ids', type=str, help='path of output file for read ids')
args = parser.parse_args()


# File opening
# *********
#ip positions
ipdb = vx.read_csv(args.ip_file, sep="\t")
#reads
reads = vx.open(args.bed_gtf_filtered)


# filter ip position with exon only present among the reads (1)
ex_tofilter = (reads[['exon_id']].to_pandas_df().astype({'exon_id':'str'})).exon_id.unique().tolist()
ipdb ['exon_id'] = ipdb.exon_id.astype("str")
ipdb = ipdb[ipdb.exon_id.isin(ex_tofilter)]

if len(ipdb) != 0:
    # Annotate reads with internal priming coordinates limits(2)
    ip_plus = ipdb[ipdb.strand=="+"].groupby("transcript_name").agg({'start_ip':'max'})
    ip_plus.rename("start_ip","threshold_coord_ip")
    ip_neg = ipdb[ipdb.strand=="-"].groupby("transcript_name").agg({'end_ip':'min'})
    ip_neg.rename("end_ip","threshold_coord_ip")
    ipdb = vx.concat([ip_plus, ip_neg])


    # Delete all the reads at an upstream position of an internal priming interval (3)
    reads = reads.join(ipdb, on="transcript_name", how="left")
    reads['fidt'] = reads.frag_id + reads.transcript_name
    to_delete1 = reads[((reads.strand == "+") & (reads.end_rd<reads.threshold_coord_ip))][["fidt"]].to_pandas_df().fidt.tolist()
    to_delete2 = reads[((reads.strand == "-") & (reads.start_rd>reads.threshold_coord_ip))][["fidt"]].to_pandas_df().fidt.tolist()
    reads = reads[~reads.fidt.isin(to_delete1)]
    reads = reads[~reads.fidt.isin(to_delete2)]
    reads = reads.drop(["fidt","threshold_coord_ip"])



def drop_duplicates(df, columns=None):
    """Return a :class:`DataFrame` object with no duplicates in the given columns.
    .. warning:: The resulting dataframe will be in memory, use with caution.
    :param columns: Column or list of column to remove duplicates by, default to all columns.
    :return: :class:`DataFrame` object with duplicates filtered away.
    """
    if columns is None:
        columns = df.get_column_names()
    if type(columns) is str:
        columns = [columns]
    return df.groupby(columns, agg={'__hidden_count': vx.agg.count()}).drop('__hidden_count')


#get final reads ids to save for the bam filtering (4)
# read_ids = drop_duplicates(reads, columns="read_id")
read_ids = reads[['read_id','transcript_name']]


#select final columns to keep
reads = reads[["seqname", "start_rd", "end_rd", "strand", "bc", "umi","start","end","gene_name","gene_id","transcript_name","transcript_id","features","exon_number","rel_start","rel_end","salmon_rlp","rel_start_rd","rel_end_rd","dist_END"]]

# Filter reads associated to genes with unique transcripts (5)
reads_unique = reads[["gene_name","transcript_name"]].to_pandas_df().drop_duplicates().groupby("gene_name").count().reset_index(drop=False)
reads_unique = reads_unique[reads_unique.transcript_name == 1]
reads_unique = reads[reads.gene_name.isin(reads_unique.gene_name)]


#write files
reads.export_csv(args.output_path, sep="\t", index=False, doublequote=False, header=False)
reads_unique.export_csv(args.output_bed_unique, sep="\t", doublequote=False, index=False, header=False)
read_ids.export_csv(args.read_ids, sep="\t", index=False, doublequote=False, header=False)

