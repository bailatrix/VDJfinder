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

# Motif description: Trp and Ser-Ser separated by exactly 8 aa
# Parameter notes: Nonamer matches >=5, heptamer matches >= 5 and no restriction on sum finds all J genes with no false positives
default = {
    'IGH': {
        'motif': 'W[A-Z]{8}SS',
        'heptamer_consensus': 'cactgtg',
        'nonamer_consensus': 'ggtttttgt',
        'gap': 23,
        'heptamer_match_min': 5,
        'nonamer_match_min': 5,
        'total_match_min': 0
        },
    
    'IGK': {
        'motif': '[FW][AG][A-Z]GT[KR][LV][DE]IK',
        'heptamer_consensus': 'cactgtg',
        'nonamer_consensus': 'ggtttttgt',
        'gap': 23,
        'heptamer_match_min': 5,
        'nonamer_match_min': 5,
        'total_match_min': 0
        },
    
    'IGL': {
        'motif': 'FG[A-Z]GT[EKQ][LV][A-Z]{2}L',
        'heptamer_consensus': 'cactgtg',
        'nonamer_consensus': 'ggtttttgt',
        'gap': 12,
        'heptamer_match_min': 5,
        'nonamer_match_min': 5,
        'total_match_min': 0
        },
    
    # NOTE: TRA and TRB loci models have not yet been implemented
    # The following two dictionaries are temporary placeholders
    'TRA': {
        'motif': 'W[A-Z]{8}SS',
        'heptamer_consensus': 'cactgtg',
        'nonamer_consensus': 'ggtttttgt',
        'gap': 0,
        'heptamer_match_min': 0,
        'nonamer_match_min': 0,
        'total_match_min': 0
        },
    
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

# Default parameters should be used for gene searches in external programs
# Set warnings for screening user-input in the case that default parameters are not used 
warnings = {
    'hept_low': "Minimum match of at least 5 finds all J genes with no false positives.",
    'hept_high': "Increasing the minimum heptamer match can lead to loss of real genes.",
    'nona_low': "Minimum match of at least 5 finds all J genes with no false positives.",
    'nona_high': "Increasing the minimum nonamer match can lead to loss of real genes.",
    'total_match': "No minimum needed for total match."
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

class J_gene ():
    
    def __init__ ( gene, locus, heptamer, nonamer, rules_dict ):
        self.gene = gene
        self.locus = locus
        self.heptamer = heptamer
        self.nonamer = nonamer
        self.rules = rules_dict
        self.matches = {}
    
    def set_matches ( self, hept_match, nona_match, total_match ):
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
def j_search( locus, jgene_file, last_v_nt=False, pseudogenes=False, overwrite=False ):    
    pseudogene_list = []
    
    # Prepare output file and build local reference database
    j_file, mode = prep_IO.prep_output( jgene_file, force=overwrite )
    jout = open( j_file, mode )
    jseq_dict, jtype_dict = prep_IO.prep_database( jgene_file )
    jout.write("> Allele | Gene Length | Start Nucleotide | End Nucleotide | Heptamer | Match | Nonamer | Match | Spacer | Total Match | Notes |\n\n\n")
    
    if pseudogenes:
        pseudo_file, mode = prep_IO.prep_pseudo_file( jgene_file, force=overwrite )
    
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
    
        # Search for W - 8 residues - SS
        # Then search backwards to conserved flanking heptamer
    
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