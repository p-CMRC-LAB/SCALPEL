

#import libs
import argparse
import pandas as pd


#Parse args
#**********
parser = argparse.ArgumentParser(description='Annotate Internal priming position')
parser.add_argument('ips', type=str, help='path of ip file')
parser.add_argument('exons', type=str, help='path of exon file')
parser.add_argument('threshold_dist', type=int, help='distance of ip position from isoform end')
parser.add_argument('output', type=str, help='path of output file')
args = parser.parse_args()

THRESHOLD_DIST = args.threshold_dist

#Read file
ipdb = pd.read_csv(args.ips, sep="\t", header=None, names=['seqname','start_ip','end_ip','strand','exon_id'])
exons = pd.read_csv(args.exons, sep="\t")

#remove NaN occurences and annotate ip occurences
ipdb = ipdb[~ipdb.exon_id.isnull()]
#explode exon_column
ipdb.exon_id = ipdb.exon_id.str.split(";")
ipdb = ipdb.explode("exon_id").drop_duplicates()

#Annotate ip positions with gene and transcript names and coords
ipdb = pd.merge(ipdb, exons[['seqname','start','end','strand','gene_name','transcript_name','exon_id','exon_ind']], on=['seqname','strand','exon_id'], how='left').dropna()

#discard ip position close to an isoform end based on a threshold
ipdb["id"] = ipdb.index
to_discard_neg = ipdb[(ipdb.exon_ind==0) & (ipdb.strand=="-") & (ipdb.start_ip < ipdb.start + THRESHOLD_DIST)].id
to_discard_pos = ipdb[(ipdb.exon_ind==0) & (ipdb.strand=="+") & (ipdb.end_ip < ipdb.end - THRESHOLD_DIST)].id
ipdb = ipdb[~ipdb.id.isin(to_discard_neg)]
ipdb = ipdb[~ipdb.id.isin(to_discard_pos)]
ipdb = ipdb.drop("id", axis=1)

#write ipfile
ipdb.to_csv(args.output, sep="\t", index=False)






