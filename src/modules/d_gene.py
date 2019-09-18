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

# Motif description: Due to short length and lack of highly conserved regions, motif needs to include heptamers
# Motif description: Look for 10-37 bases flanked by three conserved nt on each side plus additional conserved nt
default = { 
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

# Default parameters should be used for gene searches in external programs
# Set warnings for screening user-input in the case that default parameters are not used 
# Potential problems: Too low = false positives, too high = loss of real genes
warnings = {
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
    }

# Columns to organize values collected during search for output file
out_columns = {
    'allele': 'Allele',
    'gene_len': 'Gene Length',
    'st_nt': 'Start Nucleotide',
    'end_nt': 'End Nucleotide',
    'up_hept': 'Upstream Heptamer',
    'up_hept_match': 'Upsteam Heptamer Match',
    'up_nona': 'Upstream Nonamer',
    'up_nona_match': 'Upstream Nonamer Match',
    'up_spacer': 'Upstream Spacer',
    'up_total_match': 'Upstream Total Match',
    'down_hept': 'Downstream Heptamer',
    'down_hept_match': 'Downstream Heptamer Match',
    'down_nona': 'Downstream Nonamer',
    'down_nona_match': 'Downstream Nonamer Match',
    'down_spacer': 'Downstream Spacer',
    'down_total_match': 'Downstream Total Match',
    'total_match': 'Total Match',
    'notes': 'Notes'
   }

class IGHD_gene ():
    default_rules = default
    
    def __init__ ( gene, heptamer, nonamer ):
        self.gene = gene
        self.heptamer = heptamer
        self.nonamer = nonamer
        self.matches = {}
        self.rules = {}
    
    def set_matches ( self, up_hept_match, up_nona_match, up_total_match,
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
        
    def set_rules ( self, motif,
                   up_hept_cons, up_nona_cons, down_hept_cons,  down_nona_cons,
                   up_hept_match_min, up_nona_match_min, down_hept_match_min, down_nona_match_min, total_match_min ):
        self.rules = {motif,
                        up_hept_cons
                        up_nona_cons,
                        down_hept_cons, 
                        down_nona_cons, 
                        up_hept_match_min,
                        up_nona_match_min,
                        down_hept_match_min,
                        down_nona_match_min,
                        total_match_min
                     }
        
    def get_matches ( self ):
        return self.matches
    
    def get_rules ( self ):
        return self.rules

# Search for d genes in given dgene_file based on given locus file
def d_search( locus, dgene_file, last_v_nt=False, pseudogenes=False, overwrite=False ):    
    pseudogene_list = []
    
    # Prepare output file and build local reference database
    d_file, mode = prep_IO.prep_output( dgene_file, force=overwrite )
    dout = open( d_file, mode )
    dseq_dict, dtype_dict = prep_IO.prep_database( dgene_file )
    dout.write(f"> {out_columns} |\n\n\n")
    
    if pseudogenes:
        pseudo_file, mode = prep_IO.prep_pseudo_file( dgene_file, force=overwrite )
    
    # Read fasta sequence
    for nt in SeqIO.parse(locus, "fasta"):
        nts = str(nt.seq).lower()
    
        # Get nt sequences in three reading frames
        nt_len = len(nt.seq)
        nt_frame = ['','','']
        nt_frame[0] = nt[0 : nt_len - (nt_len + 0)%3]
        nt_frame[1] = nt[1 : nt_len - (nt_len + 2)%3]
        nt_frame[2] = nt[2 : nt_len - (nt_len + 1)%3]
    
        # Get aa translations in three reading frames
        aa_frame = ['','','']
        aa_frame[0] = str(nt_frame[0].seq.translate())
        aa_frame[1] = str(nt_frame[1].seq.translate())
        aa_frame[2] = str(nt_frame[2].seq.translate())
    
        # Read fasta sequence
        for dseq in finditer( ighd_motif, nts ):
            sta_nt = dseq.span()[0]
            end_nt = dseq.span()[1]
            gene = nts[sta_nt+7:end_nt-7]
            upstream_heptamer = nts[sta_nt:sta_nt+7]
            upstream_spacer = nts[sta_nt-12:sta_nt]
            upstream_nonamer = nts[sta_nt-12-9:sta_nt-12]
            downstream_heptamer = nts[end_nt-7:end_nt]
            downstream_spacer = nts[end_nt:end_nt+12]
            downstream_nonamer = nts[end_nt+12:end_nt+12+9]
            gene_len = end_nt - sta_nt - 14
    
            upstream_heptamer_match = sum(c1 == c2 for c1, c2 in zip(upstream_heptamer, ighd_upstream_heptamer_consensus))
            upstream_nonamer_match = sum(c1 == c2 for c1, c2 in zip(upstream_nonamer, ighd_upstream_nonamer_consensus))
            downstream_heptamer_match = sum(c1 == c2 for c1, c2 in zip(downstream_heptamer,ighd_downstream_heptamer_consensus))
            downstream_nonamer_match = sum(c1 == c2 for c1, c2 in zip(downstream_nonamer, ighd_downstream_nonamer_consensus))
    
            upstream_total_match = upstream_heptamer_match + upstream_nonamer_match
            downstream_total_match = downstream_heptamer_match + downstream_nonamer_match
            total_match = upstream_total_match + downstream_total_match
    
            # Look for exact matches between discovered gene and known alleles
            # Allow for the possibility that multiple alleles have same sequence
            # Default assumption is that gene is not in database
    
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
    
            if upstream_heptamer_match >= ighd_min_upstream_heptamer_match              \
                and upstream_nonamer_match >= ighd_min_upstream_nonamer_match       \
                and downstream_heptamer_match >= ighd_min_downstream_heptamer_match \
                and downstream_nonamer_match >= ighd_min_downstream_nonamer_match   \
                and total_match >= ighd_min_total_match:
            
                dout.write(f"> {allele} | {gene_len} | {sta_nt} | {end_nt} | {upstream_heptamer} | {upstream_heptamer_match} | {upstream_nonamer} | {upstream_nonamer_match} | {upstream_spacer} | {upstream_total_match} | {downstream_heptamer} | {downstream_heptamer_match} | {downstream_nonamer} | {downstream_nonamer_match} | {downstream_spacer} | {downstream_total_match} | {total_match} | {notes} |\n")
                dout.write(gene)
                dout.write("\n\n")  
            else:
                pseudogene_list.append(gene)
    
    dout.close()
    
    if pseudogenes:
        with open( pseudo_file, mode ) as pickle_file:
            pickle.dump( pseudogene_list, pickle_file )
        pickle_file.close()
    
    if last_v_nt:
        return last_v_nt