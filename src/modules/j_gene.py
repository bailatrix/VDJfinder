# system dependencies
from argparse import ArgumentParser
from os.path import basename
from Bio.Seq import Seq
from Bio import SeqIO
import pickle

# program dependencies
from re import finditer
from re import match

# homemade programs
from modules import prep_IO

# Note from the Programmer:
# These values have been rigorously tested and proven to identify and return the maximum correct immunoglobulin genes. Any alterations to the rules below risks compromising the accuracy of the program. 
# If necessary, custom search criteria can be run temporarily through the main search() method. It is not necessary to change this code.
default = {
        
    # *** DO NOT CHANGE THESE VALUES. ***
    'IGH': {
        'motif': 'W[A-Z]{8}SS',
        'hept_cons': 'cactgtg',
        'nona_cons': 'ggtttttgt',
        'gap': 23,
        'hept_match_min': 5,
        'nona_match_min': 5,
        'total_match_min': 0
        },
        
    # *** DO NOT CHANGE THESE VALUES. ***
    'IGK': {
        'motif': '[FW][AG][A-Z]GT[KR][LV][DE]IK',
        'hept_cons': 'cactgtg',
        'nona_cons': 'ggtttttgt',
        'gap': 23,
        'hept_match_min': 5,
        'nona_match_min': 5,
        'total_match_min': 0
        },
        
    # *** DO NOT CHANGE THESE VALUES. ***
    'IGL': {
        'motif': 'FG[A-Z]GT[EKQ][LV][A-Z]{2}L',
        'hept_cons': 'cactgtg',
        'nona_cons': 'ggtttttgt',
        'gap': 12,
        'hept_match_min': 5,
        'nona_match_min': 5,
        'total_match_min': 0
        },
        
    # *** DO NOT CHANGE THESE VALUES. ***
    'TRA': {
        'motif': 'W[A-Z]{8}SS',
        'hept_cons': 'cactgtg',
        'nona_cons': 'ggtttttgt',
        'gap': 0,
        'hept_match_min': 0,
        'nona_match_min': 0,
        'total_match_min': 0
        },
        
    # *** DO NOT CHANGE THESE VALUES. ***
    'TRB': {
        'motif': 'W[A-Z]{8}SS',
        'hept_cons': 'cactgtg',
        'nona_cons': 'ggtttttgt',
        'gap': 0,
        'hept_match_min': 0,
        'nona_match_min': 0,
        'total_match_min': 0
        }
    }

# Columns to organize values collected during search for output file
out_columns = out_columns = "Allele | Gene Length | Start Nucleotide | End Nucleotide | Heptamer | Match | Nonamer | Match | Spacer | Total Match | Notes |"

class J_Gene ():
    
    def __init__ ( gene, locus, heptamer, nonamer, rules_dict, matches ):
        self.gene = gene
        self.locus = locus
        self.heptamer = heptamer
        self.nonamer = nonamer
        self.rules = rules_dict
        self.matches = set_matches( matches )

    def __str__( self ):
        return f'{gene}'
        
    def __repr__( self ):
        return f'J Gene found in a(n) {locus} .fasta file.'
    
    def __set_matches ( self, hept_match, nona_match, total_match ):
        self.matches = {
            "Heptamer Match": hept_match,
            "Nonamer Match": nona_match,
            "Total Match": total_match
                       }
        
    def get_matches ( self ):
        return self.matches
    
    def get_rules ( self ):
        return self.rules


# Search for j genes in given dgene_file based on given locus file
def j_search( locus_file, locus_type, last_v_nt, pseudogenes, file_name, force, rev_dir=False ):    
    pseudogenes_list = []
    rules = default[locus_type]
    
    # Prepare output file(s)
    j_file, j_mode = prep_IO.prep_output( locus_file, locus_type, 'J', file_name, force )
    pseudo_file, mode = prep_IO.prep_output_pseudo( locus_file, locus_type, 'J', file_name, force )
        
    # Start output file and build local reference database
    jout = open( j_file, j_mode )
    jseq_dict, jtype_dict = prep_IO.prep_database( locus_type, 'J' )
    jout.write( f'> {out_columns} \n\n' )
            
    # Read fasta sequence
    for nt in SeqIO.parse( locus_file, "fasta" ):

        if rev_dir:
            nt = nt.reverse_complement()
            
        nts = str(nt.seq).lower()
        aa_frame = prep_IO.prep_frame( nt )
        
        # Break down fasta sequence
        for frame,offset in zip( aa_frame,[0,1,2] ):
            for mcons in finditer( rules[ 'motif' ], frame ):
                sta_aa = mcons.span()[0]         # Starting aa in frame
                end_aa = mcons.span()[1]         # End aa in frame
                sta_nt = sta_aa*3 + offset       # Start nt in original seq
                end_nt = end_aa*3 + offset       # End nt in original sequence
                
                # Search upstream for heptamer
                found = False
                gap = rules[ 'gap' ]
                
                for i in range( 0,39 ):
                    pos = sta_nt - i
                    heptamer = nts[ pos-7 : pos ]
                    nonamer = nts[ pos-7-9-gap : pos-7-gap ]
                    
                    # Summarize matches by each
                    heptamer_match = sum( c1 == c2 for c1, c2 in zip( heptamer, rules[ 'hept_cons' ] ) )
                    nonamer_match = sum( c1 == c2 for c1, c2 in zip( nonamer, rules[ 'nona_cons' ] ) )
                    
                    # Calculate total matches
                    total_match = heptamer_match + nonamer_match
    
                    # Verify matches meet minimum requirements
                    if heptamer[4:7] == rules[ 'hept_cons' ][4:7]     \
                    and heptamer_match >= rules[ 'hept_match_min' ]            \
                    and nonamer_match >= rules[ 'nona_match_min' ]              \
                    and total_match >= rules[ 'total_match_min' ]:
                        sta_nt -= i
                        found = True
                        break
    
                gene     = nts[ sta_nt : end_nt ]
                gene_len = end_nt - sta_nt
    
                spacer = ''
                notes = ''
            
                # Is 'aa' meant to be the same 'aa' as in v_gene search?
                # Test J gene for stop codons
                #if '*' in aa:
                #    notes = 'Contains stop codon'
                #else:
                #    notes = ''
    
                # Look for exact matches between discovered gene and known alleles
                #     - Default assumption is that gene is not known
                # NOTE: Must allow for the possibility that multiple alleles have same sequence
                allele = ''
                for jallele, jseq in jseq_dict.items():
                    if (gene == jseq) or (gene in jseq) or (jseq in gene):
                        if allele:
                            allele += ', ' + jallele + ' ' + jtype_dict[jallele]
                        else: 
                            allele = jallele + ' ' + jtype_dict[jallele]
    
                    if gene in jseq and gene != jseq:
                        if notes:
                            notes += ', ' + 'Subset of ref seq'
                        else:
                            notes = 'Subset of ref seq'
    
                    if jseq in gene and gene != jseq:
                        if notes:
                            notes += ', ' + 'Superset of ref seq'
                        else:
                            notes = 'Superset of ref seq'
                
                if not allele:
                    allele = 'Not in J ref database'
                
                # Flag J genes that start before last V gene
                if last_v_nt and sta_nt < last_v_nt:
                    continue
                    allele += ' / located in V gene region'
    
                if found:
                    jout.write(f'> {allele} | {gene_len} | {sta_nt} | {end_nt} | {heptamer} | {heptamer_match} | {nonamer} | {nonamer_match} | {spacer} | {total_match} | {notes} |\n')
                    jout.write(gene)
                    jout.write('\n\n')
                else:
                    pseudogenes_list.append(gene)
    
    jout.close() 
    
    if pseudogenes_list:
        with open( pseudo_file, mode ) as pickle_file:
            pickle.dump( pseudogenes_list, pickle_file ) 
        pickle_file.close()
        
    if last_v_nt:
        return last_v_nt
    
def custom_j_search( locus_file, locus_type, force, last_v_nt=True ):
    return 0