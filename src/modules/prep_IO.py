from os.path import basename
from os.path import isfile
from os.path import getsize
from Bio.SeqIO import parse

def prep_pseudo_file( gene_file, force=False ):
    out_name = basename(gene_file).replace('.fasta', '_pseudogenes_found.pickle')
    out_file = f'./output/{out_name}'
    
    if force:
        return out_file, 'wb'
    elif isfile(out_file) and getsize(out_file) > 0:
        user = input(f"There is already a {out_name} file. Overwrite existing output? [y/n]\t")
        
        if user == 'y':
            print("Deleting existing file contents.")
            return out_file, 'wb'
        else:
            print("Appending runtime contents to existing file.")
            return out_file, 'ab'
    else:
        return out_file, 'wb'
    
def prep_output( gene_file, force=False ):
    out_name = basename(gene_file).replace('.fasta', '_found.fasta')    
    out_file = f'./output/{out_name}'
    
    if force:
        return out_file, 'w'
    elif isfile(out_file) and getsize(out_file) > 0:
        user = input(f"There is already a {out_name} file. Overwrite existing output? [y/n]\t")
        
        if user == 'y':
            print("Deleting existing file contents.")
            return out_file, 'w'
        else:
            print("Appending runtime contents to existing file.")
            return out_file, 'a'
    else:
        return out_file, 'w'

def prep_database( gene_file ):
    seq_dict   = {} # Dictionary of D gene sequences (key = allele, value = sequence)
    type_dict  = {} # Dictionary of D gene types (key = allele, value = gene type)

    for nt in parse(gene_file, "fasta"):
        nts    = str(nt.seq).lower()
        p = nt.description.split('|')
        allele = p[1]
        gene_type   = p[3]
        seq_dict[allele]  = nts
        type_dict[allele] = gene_type
        
    return seq_dict, type_dict