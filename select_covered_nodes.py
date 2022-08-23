# usage: python select_covered_nodes.py -d outdir -s species -o outfile [-p length_covered] [-n min_coverage]

import pandas as pd
import sys
import itertools
import re
import os
import argparse



def cigar_generator(cigar_string):
    '''
    Generates blocks: length, operation from the cigar string.
    
    Arguments:
        cigar_string (str): a sequence of of base lengths and the associated operation (CIGAR string from the alignment)
        
    Returns:
        int: number of base pairs
        str: operation [MIDNSHP=X]
    
    >>>next(cigar_generator(5=1X10=2I3=))
    (5,"=")
    '''
    cigar_string=cigar_string.lstrip('cg:Z:')
    cigar_list=re.split('(\d+[MIDNSHP=X])',cigar_string)
    for element in cigar_list:
        if element:
            yield int(element[:-1]),  element[-1]
            
def read_graph(gfa_file):
    '''
    Reads sequences from graph gfa file to dictionary.
    
    Arguments:
        gfa_file (str): path to graph file in GFAv1 format
    
    Returns:
        segments (dict): dictionary with node ids as keys and sequences as values
    '''
    segments=dict()
    with open(gfa_file, 'r') as gfa:
        for line in gfa:
            if line.startswith('S'):
                node,seq=line.split('\t')[1:3]
                segments[node]=seq.strip()
    return segments
            
def list_to_file(li,file):
    '''
    Writes list to file with one element per line.
    
    Arguments:
        li (list)
        file (str): path to file
    '''
    with open(file, 'w') as f:
        for el in li:
            f.write(el+'\n')



# arg parser
parser = argparse.ArgumentParser(description='Select nodes covered by nmin contigs and write them to file')
parser.add_argument('-a','--alignment', metavar='gaf',required=True,help='path alignment file in gaf format')
parser.add_argument('-g','--graph', metavar='gfa',required=True,help='path to variation graph in gfa format')
parser.add_argument('-o','--outfile', metavar='txt',required=True,help='path to output file')

parser.add_argument('-p','--length_covered', metavar='0.1',type=float,default=0.1,help='minimal covered node length (proportion) to increase node\'s coverage (default 0.1)')
parser.add_argument('-n','--nmin', metavar='1',default=1,type=float,help='minimal coverage (int) for a node to retain (default 1)')
args = parser.parse_args()


# read paths to data
gaf=vars(args)['alignment']
gfa_file=vars(args)['graph']
outfile=vars(args)['outfile']
cov_limit=vars(args)['length_covered']
nmin=vars(args)['nmin']

# read nodes sequences from a graph
segments=read_graph(gfa_file)

# read gaf alignment file (from GraphAligner) into pandas dataframe
df=pd.read_csv(gaf,sep='\t', header=None)
df.columns=['contig', 'contig_length', 'q_start','q_end','q_strand','path','path_length','r_start','r_end','nmatches',
                  'alignment_length', 'mapq','mismatches','ali_score','divergence','identity', 'cigar']
df=df.dropna()
# cut line to a proper format (if file was corrupted)
df['cigar']=df['cigar'].str.rstrip('0123456789')
# filter alignments by alignment score >=50
df['ali_score']=pd.to_numeric(df['ali_score'].str.replace('AS:f:',''))
df=df[df['ali_score']>=50]
# split alignment path into separate nodes list
df['nodes']=df['path'].str.lstrip('>|<').str.split('>|<')
nodes=itertools.chain(*list(df['nodes']))


# create a dictionary to store nodes coverage by contigs (node id: set of contigs aligned to it)
di={key:set() for key in segments.keys()}

# iterate for each alignment in dataframe to count node coverages
for i in range(len(df)):
    path=df.iloc[i]['nodes']
    begin=df.iloc[i]['r_start']
    node_length=len(segments[path[0]])-begin # first node length
    node_cov=0 # coverage of a node by alignment (bp)
    ref_track=0 # current position on the reference
    cigar_string=df.iloc[i]['cigar']
    cigar=cigar_generator(cigar_string)
    
    # count number of mathes on the node according to cigar
    for node in path:
        node_length=len(segments[node])
        while ref_track < node_length:
            #print(nbp,operation,node, '; ref_track', ref_track, 'node_length',node_length)
            try:
                nbp,operation=next(cigar)
            except StopIteration:
                ref_track=node_length
                continue
            if operation == '=':
                ref_track+=nbp
                node_cov+=nbp
            elif operation == 'I' or operation=='X':
                ref_track+=nbp
            elif operation == 'D':
                pass
            else:
                print('Unknown operation', operation, 'in CIGAR')
                break
                
        # count coverage and decide if we keep the node
        count_cov=node_cov
        if node_cov > node_length:
            # excess in cigar=ref_track-node_length
            count_cov=node_cov - (ref_track - node_length)
        cov=count_cov/node_length
        if cov >= cov_limit:
            di[node].add(df.iloc[i]['contig'])
            #selected_nodes.append(node)

        # update coverage and tracker
        if node_cov > node_length:
            node_cov=ref_track-node_length
        ref_track=ref_track-node_length

filt_nodes=[]
for key,value in di.items():
    if len(value)>=nmin or len(segments[key])<50:
        filt_nodes.append(key)

list_to_file(filt_nodes,outfile)
