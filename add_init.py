# usage: contig_ids | python add_init.py -i init_bin_path -s sample_path | contig_ids

import pandas as pd
import statistics as st
import sys
import argparse



def headers_from_fasta(fasta):
    headers=list()
    with open(fasta) as f:
        for line in f:
            if line.startswith('>'):
                headers.append(line.strip()[1:])
    return headers

def fasta_to_df(fasta):
    headers=headers_from_fasta(fasta)
    df=pd.DataFrame(headers).iloc[:,0].str.split(' ',expand=True)
    df.columns=['contig_name','flag','coverage','length']
    df['coverage']=pd.to_numeric(df['coverage'].str.lstrip('multi='))
    return df



# arg parser
parser = argparse.ArgumentParser(description='Add contig ids from the initial bin to ids from stdin')
parser.add_argument('-i','--init_bin', metavar='path1',required=True,help='path to initial bin')
parser.add_argument('-s','--sample', metavar='path2',required=True,help='path to sample (unzipped fasta megahit-metabat output)')
args = parser.parse_args()

# read data
sample=fasta_to_df(vars(args)['sample'])
init_bin=headers_from_fasta(vars(args)['init_bin'])
ref_bin=[_id.strip().strip('"') for _id in sys.stdin.readlines()]
init_bin_df=sample[sample['contig_name'].isin(init_bin)]
ref_bin_df=sample[sample['contig_name'].isin(ref_bin)]

# filtering could be here

# add contigs from initial bin and write to stdout

filtered_names=set(ref_bin)
if len(filtered_names) != 0:
    for _id in init_bin_df['contig_name']:
        filtered_names.add(_id)
else:
    print("ERROR: There are no contigs in the refined bin. Terminating execution",file=sys.stderr)
    sys.exit(1)

for _id in filtered_names:
    print(_id)
