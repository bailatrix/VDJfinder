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
        'motif': 'C[A-Z]{8,17}W[A-Z]{38,47}[IVLFCMA][A-Z]{12,14}C',
        'hept_cons': 'cacagtg',
        'nona_cons': 'acaaaaacc',
        'gap': 23,
        'nt_before_cys1': 63,
        'hept_match_min': 5,
        'nona_match_min': 6,
        'total_match_min': 13
        },
    
    # *** DO NOT CHANGE THESE VALUES. ***
    'IGK': {
        'motif': 'C[A-Z]{8,17}W[A-Z]{38,47}[IVLFCMA][A-Z]{12,14}C',
        'hept_cons': 'cacagtg',
        'nona_cons': 'acaaaaacc',
        'gap': 12,
        'nt_before_cys1': 66,
        'hept_match_min': 5,
        'nona_match_min': 6,
        'total_match_min': 13
        },
    
    # *** DO NOT CHANGE THESE VALUES. ***
    'IGL': {
        'motif': 'C[A-Z]{8,17}W[A-Z]{38,47}[IVLFCMA][A-Z]{12,14}C',
        'hept_cons': 'cacagtg',
        'nona_cons': 'acaaaaacc',
        'gap': 23,
        'nt_before_cys1': 63,
        'hept_match_min': 5,
        'nona_match_min': 5,
        'total_match_min': 12
        },
    
    # *** DO NOT CHANGE THESE VALUES. ***
    'TRA': {
        'motif': 'C[A-Z]{8,17}W[A-Z]{38,47}[IVLFCMA][A-Z]{12,14}C',
        'hept_cons': 'cacagtg',
        'nona_cons': 'acaaaaacc',
        'gap': 0,
        'nt_before_cys1': 0,
        'hept_match_min': 0,
        'nona_match_min': 0,
        'total_match_min': 0
        },
    
    # *** DO NOT CHANGE THESE VALUES. ***
    'TRB': {
        'motif': 'C[A-Z]{8,17}W[A-Z]{38,47}[IVLFCMA][A-Z]{12,14}C',
        'hept_cons': 'cacagtg',
        'nona_cons': 'acaaaaacc',
        'gap': 0,
        'nt_before_cys1': 0,
        'hept_match_min': 0,
        'nona_match_min': 0,
        'total_match_min': 0
        }
    }

# columns to organize values collected during search for output file
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


class V_Gene ():   
    def __init__ ( gene, locus, heptamer, nonamer, rules_dict ):
        self.gene = gene
        self.locus = locus
        self.heptamer = heptamer
        self.nonamer = nonamer
        self.rules = rules_dict
        self.custom_rules = ( rules_dict == default[locus] )        # track if gene found w/ custom rules
        self.matches = {}
        
    def __str__( self ):
        return f'{gene}'
        
    def __repr__( self ):
        return f'V Gene found in a(n) {locus} .fasta file.'
    
    def __len__( self ):
        return len(self.gene)
    
    def __set_matches ( self, hept_match, nona_match, total_match ):
        self.matches = { "Heptamer Match": hept_match,
                       "Nonamer Match": nona_match,
                       "Total Match": total_match
                       }
        
    def get_matches ( self ):
        return self.matches
    
    def get_rules ( self ):
        return self.rules

# Search for v genes in given vgene_file based on given locus file
def v_search( locus_file, locus_type, last_v_nt=True, rev_dir=False ):
    pseudogene_list = []
    
    # Prepare output file and build local reference database
    v_file, mode = prep_IO.prep_output( locus_file )
    vout = open( v_file, mode )
    vseq_dict, vtype_dict = prep_IO.prep_database( locus_type, 'V' )
    vout.write("> {out_columns} \n\n\n")
    
    pseudo_file, mode = prep_IO.prep_output( locus_file, pseudogenes=True )
    
    # Read fasta sequence
    for nt in SeqIO.parse(locus_file, "fasta"):
        
        if rev_dir:
            nt = nt.reverse_complement()
            
        nts = str(nt.seq).lower()
        aa_frame = prep_IO.prep_frame( nt )
        
        for frame,offset in zip( aa_frame, [0,1,2] ):
            for mcons in finditer( ighv_motif, frame ):
                sta_aa = mcons.span()[0]         # Starting aa in frame
                end_aa = mcons.span()[1]         # End aa in frame
                sta_nt = sta_aa*3 + offset       # Start nt in original seq
                end_nt = end_aa*3 + offset       # End nt in original sequence
                
                # Search for 3' flanking sequence
                found     = False
                heptamer  = ''
                spacer    = ''
                nonamer   = ''
                post_cys2 = 0
                flank = nts[end_nt:end_nt+60]
                for i in range(0,20):
                    heptamer  = flank[i:i+7]
                    spacer  = flank[i+7:i+7+23]
                    nonamer = flank[i+7+23:i+7+23+9]
    
                    heptamer_match = sum(c1 == c2 for c1, c2 in zip(heptamer, ighv_heptamer_consensus))
                    nonamer_match = sum(c1 == c2 for c1, c2 in zip(nonamer, ighv_nonamer_consensus))
                    total_match = heptamer_match + nonamer_match

                    if heptamer[0:3] == ighv_heptamer_consensus[0:3]      \
                        and heptamer_match >= ighv_min_heptamer_match \
                        and nonamer_match >= ighv_min_nonamer_match   \
                        and total_match >= ighv_min_total_match:
                        end_nt += i
                        post_cys2 = i
                        found = True
                        break
    
                # Set start of gene to 63 bases before first conserved Cys in IGH
                # (Use 66 bases for IGK and IGL)
                sta_nt  -= ighv_nt_before_cys1
                gene     = nts[sta_nt:end_nt]
                aa       = frame[sta_aa - int(ighv_nt_before_cys1/3):end_aa + int(post_cys2/3)]
                gene_len = end_nt - sta_nt
                
                # Test V gene for stop codons
                if '*' in aa:
                    notes = 'Contains stop codon'
                else:
                    notes = ''
    
                # Look for exact matches between discovered gene and known alleles
                # Allow for the possibility that multiple alleles have same sequence
                # Default assumption is that gene is not in database
    
                allele = ''
                for vallele, vseq in vseq_dict.items():
                    if (gene == vseq) or (gene in vseq) or (vseq in gene):
                        if allele:
                            allele += ', ' + vallele + ' ' + vtype_dict[vallele]
                        else: 
                            allele = vallele + ' ' + vtype_dict[vallele]
    
                    if gene in vseq and gene != vseq:
                        if notes:
                            notes += ', ' + 'Subset of ref seq'
                        else:
                            notes = 'Subset of ref seq'
    
                    if vseq in gene and gene != vseq:
                        if notes:
                            notes += ', ' + 'Superset of ref seq'
                        else:
                            notes = 'Superset of ref seq'
                
                if not allele:
                    allele = 'Not in V ref database'
                
                if found:                
                    vout.write(f"> {allele} | {gene_len} | {sta_nt} | {end_nt} | {heptamer} | {heptamer_match} | {nonamer} | {nonamer_match} | {spacer} | {total_match} | {notes} |\n")
                    vout.write(gene)
                    vout.write("\n\n")
                    
                    last_v_nt = end_nt # Keep track of where last V gene found
                else:
                    pseudogene_list.append(gene)
    
    vout.close()
    
    if pseudogenes_list:
        with open( pseudo_file, mode ) as pickle_file:
            pickle.dump( pseudogene_list, pickle_file )
        pickle_file.close()
        
    if last_v_nt:
        return last_v_nt
    
def custom_v_search( locus_file, locus_type, last_v_nt=True ):
    return 0