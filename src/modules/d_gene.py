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
        'motif': '[acgt]ac[acgt]gtg[actg]{10,37}cac[acgt]g[actg]{2}',
        'up_hept_cons': 'cactgtg',
        'up_nona_cons': 'tgtttttgg',       #RSS question of whether last base is g or t
        'down_hept_cons': 'cacagtg',
        'down_nona_cons': 'acaaaaacc',
        'up_hept_match_min': 5,
        'up_nona_match_min': 4,
        'down_hept_match_min': 5,
        'down_nona_match_min': 4,
        'total_match_min': 22
        },
        
    # *** DO NOT CHANGE THESE VALUES. ***
    'TRB': {
        'motif': '[acgt]ac[acgt]gtg[actg]{10,37}cac[acgt]g[actg]{2}',
        'up_hept_cons': 'cactgtg',
        'up_nona_cons': 'tgtttttgg',       #RSS question of whether last base is g or t
        'down_hept_cons': 'cacagtg',
        'down_nona_cons': 'acaaaaacc',
        'up_hept_match_min': 0,
        'up_nona_match_min': 0,
        'down_hept_match_min': 0,
        'down_nona_match_min': 0,
        'total_match_min': 0
        }
    }

# Columns to organize values collected during search for output file
out_columns = out_columns = "Allele | Gene Len | Start Nt | End Nt | Upstream Hept | Match | Upstream Nona | Match | Upstream Spacer | Upstream Total Match | Downstream Hept | Match | Downstream Nona | Match | Downstream Spacer | Downstream Total Match | Total Match | Notes |"

class D_Gene ():
    
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
        return f'D Gene found in a(n) {locus} .fasta file.'
    
    def __set_matches ( self, up_hept_match, up_nona_match, up_total_match,
                     down_hept_match, down_nona_match, down_total_match, total_match ):
        self.matches = {
            'Upstream Heptamer Match': up_hept_match,
            'Upstream Nonamer Match': up_nona_match,
            'Upstream Total Match': up_total_match,
            'Downstream Heptamer Match': down_hept_match,
            'Downstream Nonamer Match': down_nona_match,
            'Downstream Total Match': down_total_match,
            'Total Match': total_match,
                       }
        
    def get_matches ( self ):
        return self.matches
    
    def get_rules ( self ):
        return self.rules

# Search for d genes in given dgene_file based on given locus file
def d_search( locus_file, locus_type, last_v_nt, pseudogenes, file_name, force, rev_dir=False ):    
    pseudogenes_list = []
    rules = default[locus_type]
    
    # Prepare output file(s)
    d_file, d_mode = prep_IO.prep_output( locus_file, locus_type, 'D', file_name, force )
    pseudo_file, mode = prep_IO.prep_output_pseudo( locus_file, locus_type, 'D', file_name, force )
    
    # Start output file and build local reference database
    dout = open( d_file, d_mode )
    dseq_dict, dtype_dict = prep_IO.prep_database( locus_type, 'D' )
    dout.write( f'> {out_columns} \n\n' )
        
    # Read fasta sequence
    for nt in SeqIO.parse( locus_file, "fasta" ):
        
        if rev_dir:
            nt = nt.reverse_complement()
            
        nts = str(nt.seq).lower()
        aa_frame = prep_IO.prep_frame( nt )
    
        # Break down fasta sequence
        for dseq in finditer( rules[ 'motif' ], nts ):
            sta_nt              = dseq.span()[0]
            end_nt              = dseq.span()[1]
            gene                = nts[ sta_nt+7 : end_nt-7] 
            upstream_heptamer   = nts[ sta_nt : sta_nt+7 ]
            upstream_spacer     = nts[ sta_nt-12 : sta_nt ]
            upstream_nonamer    = nts[ sta_nt-12-9 : sta_nt-12 ]
            downstream_heptamer = nts[ end_nt-7 : end_nt ]
            downstream_spacer   = nts[ end_nt : end_nt+12 ]
            downstream_nonamer  = nts[ end_nt+12 : end_nt+12+9 ]
            gene_len            = end_nt - sta_nt - 14
    
            # Summarize matches by each
            upstream_heptamer_match   = sum(c1 == c2 for c1, c2 in zip(upstream_heptamer, rules[ 'up_hept_cons' ] ) )
            upstream_nonamer_match    = sum(c1 == c2 for c1, c2 in zip(upstream_nonamer, rules[ 'up_nona_cons' ] ) )
            downstream_heptamer_match = sum(c1 == c2 for c1, c2 in zip(downstream_heptamer, rules[ 'down_hept_cons' ] ) )
            downstream_nonamer_match  = sum(c1 == c2 for c1, c2 in zip(downstream_nonamer, rules[ 'down_nona_cons' ] ) )
            
            # Calculate total matches
            upstream_total_match   = upstream_heptamer_match + upstream_nonamer_match
            downstream_total_match = downstream_heptamer_match + downstream_nonamer_match
            total_match            = upstream_total_match + downstream_total_match
    
            # Look for exact matches between discovered gene and known alleles
            #     - Default assumption is that gene is not known
            # NOTE: Must allow for the possibility that multiple alleles have same sequence
    
            notes = ''
            allele = 'Not in D ref db'
            for dallele, dseq in dseq_dict.items():
                if gene == dseq:
                    if allele == 'Not in D ref db':
                        allele = dallele + ' ' + dtype_dict[dallele]
                    else:
                        allele = allele + ', ' + dallele + ' ' + dtype_dict[dallele]
    
            # Flag D genes that start before last V gene
            if last_v_nt and sta_nt < last_v_nt:
                continue
                
            # Verify matches meet minimum requirements
            if upstream_heptamer_match >= rules[ 'up_hept_match_min' ]              \
                and upstream_nonamer_match >= rules[ 'up_nona_match_min' ]          \
                and downstream_heptamer_match >= rules[ 'down_hept_match_min' ]     \
                and downstream_nonamer_match >= rules[ 'down_nona_match_min' ]      \
                and total_match >= rules[ 'total_match_min' ]:
            
                dout.write(f'> {allele} | {gene_len} | {sta_nt} | {end_nt} | {upstream_heptamer} | {upstream_heptamer_match} | {upstream_nonamer} | {upstream_nonamer_match} | {upstream_spacer} | {upstream_total_match} | {downstream_heptamer} | {downstream_heptamer_match} | {downstream_nonamer} | {downstream_nonamer_match} | {downstream_spacer} | {downstream_total_match} | {total_match} | {notes} |\n')
                dout.write(gene)
                dout.write('\n\n')  
            else:
                pseudogenes_list.append(gene)
    
    dout.close()
    
    if pseudogenes_list:
        with open( pseudo_file, mode ) as pickle_file:
            pickle.dump( pseudogenes_list, pickle_file )
        pickle_file.close()
    
    if last_v_nt:
        return last_v_nt
    
def custom_d_search( locus_file, locus_type, force, last_v_nt=True ):
    return 0