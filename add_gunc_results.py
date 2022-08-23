# Usage: python add_gunc_results.py dir
# Adds a "pass GUNC" column from gunc output file to checkm2 output
#
# In the input directory there should be checkm2 and gunc folders with the corresponding output.
# Writes output to checkm2 directory: quality_report_gunc.tsv

import pandas as pd
import re
import sys


# read path of the directory with checkm2 and gunc results
dir_path=sys.argv[1]
checkm_file=f'{dir_path}/checkm2/quality_report.tsv'
GUNC_file=f'{dir_path}/gunc/GUNC.maxCSS_level.tsv'

# read checkm2 output to pandas dataframe
checkm=pd.read_csv(checkm_file,sep='\t',usecols=['Name','Completeness','Contamination'])
checkm.columns=['Bin Id','Completeness','Contamination']

# read gunc output to pandas dataframe
GUNC=pd.read_csv(GUNC_file,sep='\t',usecols=['genome','pass.GUNC'])
GUNC.columns=['Bin Id','GUNC']

# get rid of gunc refinement suffixes in bin ids (names of bins fasta files)
checkm['Bin Id']=checkm['Bin Id'].apply(lambda x: re.match('.*\.psb_metabat2\.[0-9]{5}',x).group())
GUNC['Bin Id']=GUNC['Bin Id'].apply(lambda x: re.match('.*\.psb_metabat2\.[0-9]{5}',x).group())
# merge gunc and checkm2 dataframes
checkm=pd.merge(checkm,GUNC)

# write output to checkm2 directory
checkm.to_csv(f'{dir_path}/checkm2/quality_report_gunc.tsv',sep='\t',index=False)
