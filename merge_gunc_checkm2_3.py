# Usage (from the directory with results (with gunc and checkm2 folders containing tsv output)):
# python merge_gunc_checkm.py [quality assessed directory with initial mags] [1.1 gunc refined checkm dir (optional)]
#
# Merges GUNC and checkm2 output results from initial MAGs, GUNC refined bins, and final bins into one table.

import pandas as pd
import sys
import os
import re



# usage statement: run without arguments
if len(sys.argv)<2:
  print('  Usage: \n  from the directory with results (with GUNC and CheckM2 tsv output folders):\n  python merge_gunc_checkm.py [quality assessed MAGs directory] [1.1 gunc refined checkm dir]')
  sys.exit(1)


# set variables
mags_path=sys.argv[1]
# if no path to GUNC refined (step 1.4) bins is defined, set it automatically
if len(sys.argv)==3:
  guncref_path=sys.argv[2]
else:
  guncref_path=[f.path for f in os.scandir(mags_path.replace('0_data','1.1_gunc_refined')+'/gunc/refinement/refined_genomes') if f.is_dir()][0]
new_path=os.getcwd()

checkm='checkm2'


# set paths to data: GUNC and checkm2 results from initial MAGs (init), GUNC refined bins (guncref), and final bins (new)
init_checkm_file=f'{mags_path}/{checkm}/quality_report.tsv'
new_checkm_file=f'{new_path}/{checkm}/quality_report.tsv'
init_GUNC_file=f'{mags_path}/gunc/GUNC.maxCSS_level.tsv'
new_GUNC_file=f'{new_path}/gunc/GUNC.maxCSS_level.tsv'
guncref_checkm_file=f'{guncref_path}/{checkm}/quality_report.tsv'
guncref_GUNC_file=f'{guncref_path}/gunc/GUNC.maxCSS_level.tsv'
merged_df=f'{new_path}/results.tsv'
merged_df3=f'{new_path}/results3.tsv'


# read data to pandas dataframes & set column names
init_checkm=pd.read_csv(init_checkm_file,sep='\t',usecols=['Name','Completeness','Contamination'])
new_checkm=pd.read_csv(new_checkm_file, sep='\t',usecols=['Name','Completeness','Contamination'])
guncref_checkm=(pd.read_csv(guncref_checkm_file, sep='\t',usecols=['Name','Completeness','Contamination']))

init_checkm.columns=['Bin Id','Completeness','Contamination']
new_checkm.columns=['Sample','NewCompl','NewCont']
guncref_checkm.columns=['Bin Id','GuncRefCompl','GuncRefCont']

init_GUNC=pd.read_csv(init_GUNC_file,sep='\t',usecols=['genome','pass.GUNC'])
init_GUNC.columns=['Bin Id','init_GUNC']
new_GUNC=pd.read_csv(new_GUNC_file,sep='\t',usecols=['genome','pass.GUNC'])
new_GUNC.columns=['Sample','new_GUNC']
guncref_GUNC=pd.read_csv(guncref_GUNC_file,sep='\t',usecols=['genome','pass.GUNC'])
guncref_GUNC.columns=['Bin Id','guncref_GUNC']
#
# to do: add dRep score


# merge dataframes

# transform sample names
init_checkm['Sample']=init_checkm['Bin Id'].apply(lambda x: x.split('.psa_megahit')[0])
new_checkm['Sample']=new_checkm['Sample'].apply(lambda x: x.split('-assembled')[0]).str.lstrip('mapped_')
new_GUNC['Sample']=new_GUNC['Sample'].apply(lambda x: x.split('-assembled')[0]).str.lstrip('mapped_')
# merge new checkm2 and new GUNC dataframes on sample name
new_checkm=pd.merge(new_checkm,new_GUNC,how='left')

# get rid of gunc refinement suffixes in bin ids
guncref_checkm['Bin Id']=guncref_checkm['Bin Id'].apply(lambda x: re.match('.*\.psb_metabat2\.[0-9]{5}',x).group())
guncref_GUNC['Bin Id']=guncref_GUNC['Bin Id'].apply(lambda x: re.match('.*\.psb_metabat2\.[0-9]{5}',x).group())
# merge GUNC and checkm2 results from GUNC refined bins on bin ids
guncref_checkm=pd.merge(guncref_checkm,guncref_GUNC)

# merge initial checkm2 and final checkm2&GUNC results
# how='left' -> retain multiple bins from the same sample (rely on initial bins structure)
df=pd.merge(init_checkm,new_checkm,how='left')
df=pd.merge(df,init_GUNC,how='left')

print(df.isna().any(axis=1).sum(),'bins missing')
print(df[df.isna().any(axis=1)]['Bin Id'])


# write output to files
# results.tsv - initial & final
# results3.tsv - initial & GUNC refined & final

# set column names, discard missing bins from the table, write output to file
df=df[['Bin Id','Sample','Completeness','NewCompl','Contamination','NewCont','init_GUNC','new_GUNC']]
df=df.dropna()
df.to_csv(merged_df,sep='\t',index=False)

# add results from GUNC refined bins (step 1.4) and write to file
df3=pd.merge(df,guncref_checkm,how='left')
df3=df3[['Bin Id','Sample','Completeness','GuncRefCompl','NewCompl','Contamination','GuncRefCont','NewCont','init_GUNC','guncref_GUNC','new_GUNC']]
df3=df3.dropna()
df3.to_csv(merged_df3,sep='\t',index=False)

