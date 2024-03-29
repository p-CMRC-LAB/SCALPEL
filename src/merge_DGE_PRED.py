

#import libraries
import argparse
import pandas as pd
import glob


#Argunent Parser
#**************
parser = argparse.ArgumentParser(description='Merge APA prediction with DGE Matrix')
parser.add_argument('PRED_path', metavar='pfip', type=str, help='path of cell isoform prediction file')
parser.add_argument('DGE_path', metavar='dfip', type=str, help='path of DGE file')
parser.add_argument('OUTPUT_APADGE_path', metavar='opfop', type=str, help='path of output APA_DGE file')
args = parser.parse_args()


#Open files
print('Files opening...')
predictions = pd.read_csv(args.PRED_path, names = ["bc","gene_name","transcript_name","tr_prob"], sep="\t", skiprows = 1)

print('dge...')
dge = pd.read_csv(args.DGE_path, sep="\t")
dge = dge.rename(columns = {'GENE':'gene_name'})

#melt dge table
print('melting table...')
dge = dge.melt(id_vars='gene_name')

predictions["bc_gn"] = predictions.bc + "_" + predictions.gene_name
dge["bc_gn"] = dge.variable + "_" + dge.gene_name

#filter barcodes not present in the dge
print('bc filtering..')
dge = dge[dge.bc_gn.isin(predictions.bc_gn.tolist())]
dge = dge.rename(columns = {'variable':'bc', 'value':'counts'})

#merge with predictions table
dge = predictions[['bc','gene_name','transcript_name','tr_prob']].merge(dge, on=['bc','gene_name'], how='left')

#Get relatives values on dge
dge['gene_transcript'] = dge.gene_name + '***' + dge.transcript_name
dge['relative_prob'] = dge.counts * dge.tr_prob
dge = dge[['bc','gene_transcript','relative_prob']]

#pivot of table
print('table pivoting...')
dge = dge.pivot_table(index='gene_transcript', columns='bc', values='relative_prob', fill_value=0)

#writing
print('Writing...')
dge.to_csv(args.OUTPUT_APADGE_path, sep="\t")
