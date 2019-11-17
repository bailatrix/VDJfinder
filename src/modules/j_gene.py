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
        'heptamer_consensus': 'cactgtg',
        'nonamer_consensus': 'ggtttttgt',
        'gap': 23,
        'heptamer_match_min': 5,
        'nonamer_match_min': 5,
        'total_match_min': 0
        },
        
    # *** DO NOT CHANGE THESE VALUES. ***
    'IGK': {
        'motif': '[FW][AG][A-Z]GT[KR][LV][DE]IK',
        'heptamer_consensus': 'cactgtg',
        'nonamer_consensus': 'ggtttttgt',
        'gap': 23,
        'heptamer_match_min': 5,
        'nonamer_match_min': 5,
        'total_match_min': 0
        },
        
    # *** DO NOT CHANGE THESE VALUES. ***
    'IGL': {
        'motif': 'FG[A-Z]GT[EKQ][LV][A-Z]{2}L',
        'heptamer_consensus': 'cactgtg',
        'nonamer_consensus': 'ggtttttgt',
        'gap': 12,
        'heptamer_match_min': 5,
        'nonamer_match_min': 5,
        'total_match_min': 0
        },
        
    # *** DO NOT CHANGE THESE VALUES. ***
    'TRA': {
        'motif': 'W[A-Z]{8}SS',
        'heptamer_consensus': 'cactgtg',
        'nonamer_consensus': 'ggtttttgt',
        'gap': 0,
        'heptamer_match_min': 0,
        'nonamer_match_min': 0,
        'total_match_min': 0
        },
        
    # *** DO NOT CHANGE THESE VALUES. ***
    'TRB': {
        'motif': 'W[A-Z]{8}SS',
        'heptamer_consensus': 'cactgtg',
        'nonamer_consensus': 'ggtttttgt',
        'gap': 0,
        'heptamer_match_min': 0,
        'nonamer_match_min': 0,
        'total_match_min': 0
        }
    }

# Columns to organize values collected during search for output file
out_columns = { "allele": "Allele",
               "gene_len" : "Gene Length",
               "st_nt" : "Start Nucleotide",
               "end_nt" : "End Nucleotide",
               "hept" : "Heptamer",
               "hept_match" : "Heptamer Match",
               "non" : "Nonamer",
               "non_match" : "Nonamer Match",
               "spacer" : "Spacer",
               "total_match" : "Total Match",
               "notes" : "Notes"
              }

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
def j_search( locus_file, locus_type, last_v_nt=True, rev_dir=False ):    
    pseudogene_list = []
    
    # Prepare output file and build local reference database
    j_file, mode = prep_IO.prep_output( locus_file )
    jout = open( j_file, mode )
    jseq_dict, jtype_dict = prep_IO.prep_database( locus_type, 'J' )
    jout.write(f"> {out_columns} |\n\n\n")
    
    pseudo_file, mode = prep_IO.prep_output( locus_file, pseudogenes=True )
    
    # Read fasta sequence
    for nt in SeqIO.parse(locus_file, "fasta"):

        if rev_dir:
            nt = nt.reverse_complement()
            
        nts = str(nt.seq).lower()
        aa_frame = prep_IO.prep_frame( nt )
    
        for frame,offset in zip(aa_frame,[0,1,2]):
            for mcons in finditer(ighj_motif, frame):
                sta_aa = mcons.span()[0]         # Starting aa in frame
                end_aa = mcons.span()[1]         # End aa in frame
                sta_nt = sta_aa*3 + offset       # Start nt in original seq
                end_nt = end_aa*3 + offset       # End nt in original sequence
                
                # Search upstream for heptamer
                found     = False
                for i in range(0,39):
                    heptamer = nts[sta_nt-i-7:sta_nt-i]
                    nonamer = nts[sta_nt-i-7-23-9:sta_nt-i-7-23]
                    heptamer_match = sum(c1 == c2 for c1, c2 in zip(heptamer, ighj_heptamer_consensus))
                    nonamer_match = sum(c1 == c2 for c1, c2 in zip(nonamer, ighj_nonamer_consensus))
                    total_match = heptamer_match + nonamer_match
    
                    if heptamer[4:7] == ighj_heptamer_consensus[4:7]     \
                            and heptamer_match >= ighj_min_heptamer_match \
                            and nonamer_match >= ighj_min_nonamer_match   \
                            and total_match >= ighj_min_total_match:
                        sta_nt -= i
                        found = True
                        break
    
                gene     = nts[sta_nt:end_nt]
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
                # Allow for the possibility that multiple alleles have same sequence
                # Default assumption is that gene is not in database
    
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
                    jout.write(f"> {allele} | {gene_len} | {sta_nt} | {end_nt} | {heptamer} | {heptamer_match} | {nonamer} | {nonamer_match} | {spacer} | {total_match} | {notes} |\n")
                    jout.write(gene)
                    jout.write("\n\n")
                else:
                    pseudogene_list.append(gene)
    
    jout.close() 
    
    if pseudogenes:
        with open( pseudo_file, mode ) as pickle_file:
            pickle.dump( pseudogene_list, pickle_file ) 
        pickle_file.close()
        
    if last_v_nt:
        return last_v_nt
    
def custom_j_search( locus_file, locus_type, last_v_nt=True ):
    return 0