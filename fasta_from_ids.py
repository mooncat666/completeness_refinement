import sys

def fasta_to_dict(fasta):
    # fasta file into dictionary with id as a key and sequence as a value
    di=dict()
    key='empty'
    value=''
    with open(fasta) as file:
        for line in file:
            if line.startswith('>'):
                di[key]=value
                value=''
                key=line[1:].split()[0].strip()
            else:
                value+=line.strip()
        di[key]=value
        del di['empty']
    return di

contigs=fasta_to_dict(sys.argv[1])
for _id in sys.stdin:
    _id=_id.strip()
    sys.stdout.write('>'+_id+'\n')
    sys.stdout.write(contigs[_id]+'\n')
