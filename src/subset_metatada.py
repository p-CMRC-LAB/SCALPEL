import argparse
import pandas as pd

#Argunent Parser
#-----------------
parser = argparse.ArgumentParser(description='Subset reads info columns')
parser.add_argument('bedf', metavar='reads', type=str, help='path of bed file')
parser.add_argument('output_path', metavar='Bmfip', type=str, help='path of output bed file')
args = parser.parse_args()

#import bed file
bed = pd.read_csv(args.bedf, sep="\t", header=None)

#subset column
bed = bed.iloc[:,[0,1,2,3,5,13,14]]
bed.columns = ["chr","Start","End","read_id","Strand","bc","umi"]

#add + 1 to the start coord to make it match with the 1 based GTF format
bed.Start = bed.Start + 1

bed = bed.drop_duplicates(["chr","Start","End","Strand","bc","umi"]).sort_values(["chr","Start","End"])

print(bed)

#Identify spliced reads
dict_save = {}
bed['read_id_uniq'] = bed.read_id.str.replace("/[0-9]","", regex=True)
for group in bed.groupby("read_id_uniq"):
	tab = group[1]
	if len(tab) > 1:
		dict_save[group[0]] = len(tab)
	else:
		dict_save[group[0]] = 1
A = pd.DataFrame.from_dict(dict_save, orient = 'index').reset_index()
A.columns = ["read_id_uniq","nb_splices"]

#merge nb_splice annotation into bed file
bed = bed.merge(A, on="read_id_uniq", how='left')
bed.read_id = bed.read_id_uniq
del bed["read_id_uniq"]


#write
bed.to_csv(args.output_path, sep="\t", header=False, index=False)


#NDX550569_RUO:71:HVCKLBGXK:4:13609:8701:13104
