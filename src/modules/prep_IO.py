# file dependencies
#    - used in prep_pseudo_file(), prep_output()
from os.path import basename
from os.path import isfile
from os.path import getsize

from os import getcwd

# .fasta handling dependencies
#    - used in prep_database()
from Bio.SeqIO import parse

# Locations of all reference database files
# These values should be hard-coded into the program as the location of each file is vital to the search function
all_ref_dbs = {
    'IGH': {
        'V': '/home/baileyp/projects/VDJfinder/src/data/master/IMGT/IGHV.fasta',
        'D': '/home/baileyp/projects/VDJfinder/src/data/master/IMGT/IGHD.fasta',
        'J': '/home/baileyp/projects/VDJfinder/src/data/master/IMGT/IGHJ.fasta'
        },
    
    'IGK': {
        'V': '/home/baileyp/projects/VDJfinder/src/data/master/IMGT/IGKV.fasta',
        'J': '/home/baileyp/projects/VDJfinder/src/data/master/IMGT/IGKJ.fasta'
        },
    
    'IGL': {
        'V': '/home/baileyp/projects/VDJfinder/src/data/master/IMGT/IGLV.fasta',
        'J': '/home/baileyp/projects/VDJfinder/src/data/master/IMGT/IGLJ.fasta'
        },
    
    'TRA': {
        'V': '/home/baileyp/projects/VDJfinder/src/data/master/IMGT/TRAV.fasta',
        'J': '/home/baileyp/projects/VDJfinder/src/data/master/IMGT/TRAJ.fasta'
        },
    
    'TRB': {
        'V': '/home/baileyp/projects/VDJfinder/src/data/master/IMGT/TRBV.fasta',
        'D': '/home/baileyp/projects/VDJfinder/src/data/master/IMGT/TRBD.fasta',
        'J': '/home/baileyp/projects/VDJfinder/src/data/master/IMGT/TRBJ.fasta'
        },
    }

##############################################################

# default parameters should be used for gene searches in external programs
# warnings for screening user-input in the case that default parameters are not used
warnings = {
    
    # V-genes: 1st conserved Cys to 2nd conserved Cys
    #     Max number of residues between Cys1-Trp          = 17
    #     Max number of residues between Tpr-[IVLFCMA]     = 47
    #     Max number of residues between [IVLFCMA]-Cys2    = 14
    #     Max number of missing residues in Cys1-Trp       = 9
    #     Max number of missing residues in Trp-[IVLFCMA]  = 9
    #     Max number of missing residues in [IVLFCMA]-Cys2 = 2
    'V': {
        'hept_low': "Setting the minimum heptamer match lower can increase false positives.",
        'hept_high': "Increasing the minimum heptamer match can lead to loss of real genes.",
        'nona_low': "Setting the minimum nonamer match lower can increase false positives.",
        'nona_high': "Increasing the minimum nonamer match can lead to loss of real genes.",
        'total_match': "Total matches >= 13 filters ORFs without losing real genes."
        },
    
    # Motif description: Due to short length and lack of highly conserved regions, motif needs to include heptamers
    # Motif description: Look for 10-37 bases flanked by three conserved nt on each side plus additional conserved nt
    'D': {
        'motif': 'Using an untested search motif is not recommended.',
        'up_hept_cons': '',
        'up_nona_cons': '',
        'down_hept_cons': '', 
        'down_nona_cons': '', 
        'up_hept_match_min': '',
        'up_nona_match_min': '',
        'down_hept_match_min': '',
        'down_nona_match_min': '',
        'total_match_min': ''
        },
    
    # Motif description: Trp and Ser-Ser separated by exactly 8 aa
    # Nonamer matches >=5, heptamer matches >= 5, and no restriction on sum finds all J genes with no false positives
    'J': {
        'hept_low': "Minimum match of at least 5 finds all J genes with no false positives.",
        'hept_high': "Increasing the minimum heptamer match can lead to loss of real genes.",
        'nona_low': "Minimum match of at least 5 finds all J genes with no false positives.",
        'nona_high': "Increasing the minimum nonamer match can lead to loss of real genes.",
        'total_match': "No minimum needed for total match."
        }
    } 

##############################################################
#                    SCREENING INPUTS                        #
##############################################################

# Private to public search() method
# screen_locus_file():
#    -
# Parameters:
#    -
# Return:
#    -
def screen_locus_file ( in_file ):
    try:
        print('')
    except :
        print( 'ERROR' )
        exit()

# Private to public search() method
# screen_locus_type():
#    -
# Parameters:
#    -
# Return:
#    -
def screen_locus_type ( in_locus ):
    try:
        all_ref_dbs[ in_locus ]
    except KeyError:
        print( f'ERROR:\t {in_locus} is not an available locus type.' )
        exit()
        
# Private to public search() method
# screen_():
#    -
# Parameters:
#    -
# Return:
#    -
def screen_genes ( in_locus, in_gene ):
    try:
        all_ref_dbs[ in_locus ][ in_gene ]
    except KeyError:
        print( f'ERROR:\t {in_gene} is not an available gene to search.' )
        exit()
        
# Private to public search() method
# screen_():
#    -
# Parameters:
#    -
# Return:
#    -
def screen_ ( in_ ):
    try:
        print('')
    except :
        print( 'ERROR' )
        exit()

##############################################################
#                    SCREENING OUTPUTS                       #
##############################################################

# Private to {gene}_search methods
# prep_output():
#    - preps output file to store search results
# Parameters:
#    - gene_file: file that is being searched, naming convention uses this file name as base
#    - pseudogenes: Boolean value, indicates whether to alter naming convention to include 'pseudogenes'
#    - pref_name: allows out_name to be passed as argument
#    - force: Boolean value, indicates whether to force overwrite in the event of duplicate file name
# Return: 
#    - file location+name and mode to use
def prep_output( locus_file, locus_type, gene, file_name, force ):    
    dest_dir = getcwd()+'/data/genes/'+locus_type
    
    # support for alternate naming convention for files
    if file_name: out_name = file_name    
    else:
        out_name = basename(locus_file).replace('.fasta', f'_{ locus_type+gene }_found.fasta')

    # format out_name to include file location
    out_file = f'{dest_dir}/{out_name}'
    
    # support for forced rewrite of file with duplicate name
    if force:
        return out_file, 'w'
    
    # optional user input to not rewrite duplicate file
    elif isfile(out_file) and getsize(out_file) > 0:
        user = input(f"There is already a(n) {out_name} file. Overwrite existing output? [y/n]\t")
        
        if user == 'y':
            print("Deleting existing file contents.")
            return out_file, 'w'
        else:
            print("Appending runtime contents to existing file.")
            return out_file, 'a'
    
    # out_file is unique
    else:
        return out_file, 'w'

# Private to {gene}_search methods
# prep_output():
#    - preps output file to store search results
# Parameters:
#    - gene_file: file that is being searched, naming convention uses this file name as base
#    - pseudogenes: Boolean value, indicates whether to alter naming convention to include 'pseudogenes'
#    - pref_name: allows out_name to be passed as argument
#    - force: Boolean value, indicates whether to force overwrite in the event of duplicate file name
# Return: 
#    - file location+name and mode to use
def prep_output_pseudo( locus_file, locus_type, gene, file_name, force ):    
    dest_dir = getcwd()+'/data/_unprocessed_'
    
    # support for alternate naming convention for files
    if file_name: out_name = file_name    
    else:
        out_name = basename(locus_file).replace('.fasta', f'_{ locus_type+gene }_pseudogenes_found.pickle')

    # format out_name to include file location
    out_file = f'{dest_dir}/{out_name}'
    
    # support for forced rewrite of file with duplicate name
    if force:
        return out_file, 'wb'
    
    # optional user input to not rewrite duplicate file
    elif isfile( out_file ) and getsize( out_file ) > 0:
        user = input(f"There is already a {out_name} file. Overwrite existing output? [y/n]\t")
        
        if user == 'y':
            print("Deleting existing file contents.")
            return out_file, 'wb'
        else:
            print("Appending runtime contents to existing file.")
            return out_file, 'ab'
    
    # out_file is unique
    else:
        return out_file, 'wb'
    
##############################################################
#                    PREP SEARCH DATA                        #
##############################################################
    
# Private to {gene}_search() methods
# prep_database():
#    - prepares local reference database as dictionaries
# Parameters:
#    - locus_type: 'IGH', 'IGL', 'IGK', 'TRA', or 'TRB'
#    - gene_type: 'V', 'D', or 'J'
# Return:
#    - dict of gene sequences and dict of gene types
def prep_database( locus_type, gene_type ):
    from sys import exc_info

    # attempt to build local dicts
    try:
        # screen input arguments
        locus_type = str(locus_type).upper()
        gene = str(gene_type).upper()
        
        # pull location of corresponding ref_db
        ref_db_file = all_ref_dbs[locus_type][gene_type]
        
        # prep return
        seq_dict   = {}
        type_dict  = {}
        
        # parse fasta file to build local ref_dbs as seq_dict and type_dict
        for nt in parse(ref_db_file, "fasta"):
            nts    = str(nt.seq).lower()
            p = nt.description.split('|')
            allele = p[1]
            gene_type   = p[3]
            seq_dict[allele]  = nts
            type_dict[allele] = gene_type
        
        return seq_dict, type_dict
    
    # handle incorrect number of arguments passed
    except TypeError:
        print("Invalid input. prep_database() takes exactly 2 string inputs as arguments: locus_type and gene_type.")
        raise TypeError
    
    # handle invalid arguments passed
    except ValueError:
        print("Invalid input. prep_database() takes 2 string inputs as arguments: locus_type and gene_type.")
        raise ValueError
    
    # handle missing file location for locus+gene lookup
    except KeyError:
        print( "Mismatch between locus_type and gene_type." )
        print( "Specified gene_type may not be in given locus_type, or locus_type is not included in known options." )
        raise KeyError
    
    # handle missing ref_db file for locus+gene lookup
    except FileNotFoundError:
        print( "Reference file for locus type {locus_type} and gene_type {gene_type} was not found." )
        raise FileNotFoundError
    
    # handle unknown error
    except:
        print("Unexpected error:", exc_info()[0])
        raise
        
# Private to v_gene_search() method
# prep_frame():
#    - TBA
# Parameter:
#    - nt: nucleotide sequence to operate on
# Return:
#    - amino acid frame
def prep_frame( nt ):
    nt_len = len(nt.seq)
    
    # Get nt sequences in three reading frames
    nt_frame = ['', '', '']
    nt_frame[0] = nt[0 : nt_len - (nt_len + 0)%3]
    nt_frame[1] = nt[1 : nt_len - (nt_len + 2)%3]
    nt_frame[2] = nt[2 : nt_len - (nt_len + 1)%3]
    
    # Get aa translations in three reading frames
    aa_frame = ['', '', '']
    aa_frame[0] = str(nt_frame[0].seq.translate())
    aa_frame[1] = str(nt_frame[1].seq.translate())
    aa_frame[2] = str(nt_frame[2].seq.translate())
    
    return aa_frame