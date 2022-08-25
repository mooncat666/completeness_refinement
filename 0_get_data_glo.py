# get_data.py 
# Usage: python get_data.py [species]
# species=!GTDB_tk species ('Fusobacterium nucleatum') or specI or ANI!

import pymongo
import re
import os
import sys
import subprocess
import pandas as pd



def get_path_of_bin(db, bin_id):
    '''
    Get path of a bin from global data.

    Searches for a bin by its id in Mongo DB and returns its absolute path.    
    '''
    bin_record = db.bins.find_one({'bin_id': bin_id})
    study_id = db.samples.find_one({'sample_id': bin_record['sample_id']})['study_id']
    return f"/g/scb2/bork/data/MAGs/{study_id}/psa_megahit/psb_metabat2/{bin_id}.fa.gz"

def get_path_of_reference(db, bin_id):
    path = db.samples.find_one({'sample_id': bin_id})['analysis_path']
    return f"{path}/ref_genome/{bin_id}.genomes.fa.gz"


def get_path_of_assembly(db, bin_id):
    '''
    Get path of an assembly the bin belongs to.

    Searches for a bin by its id and a sample it comes from in Mongo DB. Returns the absolute path to the assembly of this sample.
    Arguments:
        db (pymongo.database.Database): Mongo database
        bin_id (str)
        
    Returns:
        str: absolute path to assembly
    '''    
    bin_record = db.bins.find_one({'bin_id': bin_id})
    study_id = db.samples.find_one({'sample_id': bin_record['sample_id']})['study_id']
    return f"/g/scb2/bork/data/MAGs/{study_id}/psa_megahit/assemblies/{bin_record['sample_id']}-assembled.fa.gz"


def create_directory(_dir): # TO DO: implement overwrite option
    '''
    Create a directory if it doesn't exist, otherwise exit the script.
    
    Arguments: 
        _dir (str): path to directory
    '''
    if not os.path.exists(_dir):
        subprocess.run(['mkdir',_dir])
    # elif overwrite:
        # rm _dir
        # mkdir _dir
    else:
        print(f'Directory {_dir} should not exist')
        sys.exit(1)

        
def find_bins(df, species, source):
    # find bins in df with the provided specI assignment or ANI cluster name/number
    if re.match('specI.*',species):
        col='specI'
    elif re.match('ANI_\d\d_v1_\d+$',species):
        col='ANI_95'
    else:
        col='species'
    bin_ids=list(df[(df[col]==species) & (df['source']==source)]['genome_id'])
    return bin_ids


def create_links(Genus_sp, bin_ids, db):
    mags_dir=f'{os.getcwd()}/0_data/{Genus_sp}/init_MAGs'
    create_directory(mags_dir)
    
    # create links to bins in output directory
    for f in bin_ids:
        subprocess.run(['ln','-s',get_path_of_bin(db,f),mags_dir])
        
    # create links to samples assemblies in output directory
    as_dir=f'{os.getcwd()}/0_data/{Genus_sp}/samples'
    create_directory(as_dir)
    bin_ids=[re.sub('.fa(.gz)?.*','',f) for f in os.listdir(mags_dir) if os.path.isfile(os.path.join(mags_dir,f))] # get bin ids from written links

    for bin_id in bin_ids:
        subprocess.run(['ln','-s',get_path_of_assembly(db,bin_id),as_dir],stderr=subprocess.DEVNULL)

def create_links_reference(Genus_sp, bin_ids, db):
    ref_out_dir=f'{os.getcwd()}/0_data/{Genus_sp}/reference_genomes'
    subprocess.run(['mkdir','-p',ref_out_dir])
    for f in bin_ids:
        subprocess.run(['ln','-s',get_path_of_reference(db,f),ref_out_dir])

species=sys.argv[1]
try:
    db_type=str(sys.argv[2])
except IndexError:
    print('Database is not specified. Setting it to global')
    db_type='global'
    
df=pd.read_csv('/g/scb2/bork/tschmidt/global_data/data/data.genome.v02.tsv',sep='\t')
Genus_sp=sys.argv[1].replace(' ','_')
create_directory(f'{os.getcwd()}/0_data/{Genus_sp}')


client = pymongo.MongoClient('mongodb://mag_read:password_mag_read@koppa:26016')

if db_type == 'freeze13':
    source='freeze13'
    client = pymongo.MongoClient('mongodb://progenomes_read:password_progenomes_read@koppa:26017')
    db=client['progenomes']
    bin_ids=find_bins(df, species, source)
    create_links_reference(Genus_sp, bin_ids, db)
elif db_type == 'global': 
    source='glomag'
    client = pymongo.MongoClient('mongodb://mag_read:password_mag_read@koppa:26016')
    db=client['mags']
    bin_ids=find_bins(df, species, source)
    create_links(Genus_sp, bin_ids, db)
elif db_type == 'both':
    source='freeze13'
    client = pymongo.MongoClient('mongodb://progenomes_read:password_progenomes_read@koppa:26017')
    db=client['progenomes']
    bin_ids=find_bins(df, species, source)
    create_links_reference(Genus_sp, bin_ids, db)
    
    source='glomag'
    client = pymongo.MongoClient('mongodb://mag_read:password_mag_read@koppa:26016')
    db=client['mags']
    bin_ids=find_bins(df, species, source)
    create_links(Genus_sp, bin_ids, db)
else:
    print('Wrong database type. Allowed: global, freeze13, both. Exiting the script')
    sys.exit(1)



    
# get checkm results and dRep scores from global data
global_data = df
checkm_dir=f'{os.getcwd()}/0_data/{Genus_sp}/checkm'
create_directory(checkm_dir)
checkm_scores=global_data[global_data['genome_id'].isin(bin_ids)][['genome_id','completeness','contamination','drep']]
checkm_scores.columns=['Bin Id','Completeness','Contamination','dRep']
checkm_scores.to_csv(f'{checkm_dir}/stats.tsv',sep='\t',index=False)

