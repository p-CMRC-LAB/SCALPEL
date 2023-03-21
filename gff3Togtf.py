
import csv
import pandas as pd
import pyranges as pr
import argparse
import numpy as np

#Argument parser
parser = argparse.ArgumentParser(description='Conversion of gf33 into suitable GTF format for scalpel')
parser.add_argument('gff3', metavar="GFF3", type=str, help="Path of gff3 file")
parser.add_argument('output_path', metavar="GTF_name", type=str, help="Path of output gtf file")
args = parser.parse_args()

#read file
print("file opening...")
gff = pr.read_gff3(args.gff3, as_df=True)

#Modification table format
print("processing table...")
#genes annotations
gff['gene_name'] = gff.exon_id.str.replace("\\..*","", regex=True)
gff['gene_id'] = gff.gene_name
gff['gene_type'] = "protein_coding"
#transcripts annotations
gff['transcript_id'] = gff.Parent.str.replace("transcript\\:","", regex=True)
gff['transcript_name'] = gff.transcript_id
gff['transcript_type'] = "protein_coding"
gff['exon_number'] = gff.exon_id.str[-1]
gff = gff[gff.Feature=="exon"]


#Formating into GTF
print("GTF formatting...")
a1 = ['gene_id ' + '"' + str(x) + '";' for x in gff.gene_id.tolist()]
a2 = ['gene_name ' + '"' + str(x) + '";' for x in gff.gene_name.tolist()]
a3 = ['gene_type ' + '"' + str(x) + '";' for x in gff.gene_type.tolist()]
a4 = ['transcript_id ' + '"' + str(x) + '";' for x in gff.transcript_id.tolist()]
a5 = ['transcript_name ' + '"' + str(x) + '";' for x in gff.transcript_name.tolist()]
a6 = ['transcript_type ' + '"' + str(x) + '";' for x in gff.transcript_type.tolist()]
a7 = ['exon_number ' + '"' + str(x) + '";' for x in gff.exon_number.tolist()]
a8 = ['exon_id ' + '"' + str(x) + '";' for x in gff.exon_id.tolist()]
annot_list = pd.DataFrame({"c1":a1,"c2":a2,"c3":a3,"c4":a4,"c5":a5,"c6":a6,"c7":a7,"c8":a8}).values.tolist()
annot_list = [" ".join(x) for x in annot_list]
gff['annotation'] = annot_list
gff = gff[["Chromosome","Source","Feature","Start","End","Score","Strand","Frame","annotation"]]
gff = gff.replace('nan', np.nan).dropna()

#return gtf table
print("writing...")
gff.to_csv(args.output_path, sep="\t", index=False, header=False, quoting=csv.QUOTE_NONE)
