


import argparse
import pandas as pd
import vaex as vx

# Argument Parser
# -----------------
parser = argparse.ArgumentParser(description='get read_id')
parser.add_argument('ip_file', metavar='reads', type=str, help='path of ip filtered file')
parser.add_argument('output_path', metavar='Bmfip', type=str, help='path of output fragment bed file')
args = parser.parse_args()


# Open file
# ---------
reads = (vx.open(args.ip_file)).to_pandas_df()

# Write file
(reads[["read_id"]].drop_duplicates()).to_csv(args.output_path, sep="\t", index=False, header=False)
