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

# Search for j genes in given dgene_file based on given locus file
def j_search( locus, jgene_file, last_v_nt=False, pseudogenes=False, overwrite=False ):    
    pseudogene_list = []
    
    # Prepare output file and build local reference database
    j_file, mode = prep_IO.prep_output( jgene_file, force=overwrite )
    jout = open( j_file, mode )
    jseq_dict, jtype_dict = prep_IO.prep_database( jgene_file )
    
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
            for mcons in finditer('W[A-Z]{8}SS', frame):
                sta_aa = mcons.span()[0]         # Starting aa in frame
                end_aa = mcons.span()[1]         # End aa in frame
                sta_nt = sta_aa*3 + offset       # Start nt in original seq
                end_nt = end_aa*3 + offset       # End nt in original sequence
                
                # Search upstream for heptamer
                found     = False
                for i in range(0,39):
                    # upstream heptamer consensus = cactgtg
                    upstream_heptamer = nts[sta_nt-i-7:sta_nt-i]
                    upstream_heptamer_score = 0
                    if upstream_heptamer[0] == 'c': upstream_heptamer_score += 1
                    if upstream_heptamer[1] == 'a': upstream_heptamer_score += 1
                    if upstream_heptamer[2] == 'c': upstream_heptamer_score += 1
                    if upstream_heptamer[3] == 't': upstream_heptamer_score += 1
    
                    # upstream nonamer consensus = ggtttttgt
                    upstream_nonamer = nts[sta_nt-i-7-23-9:sta_nt-i-7-23]
                    upstream_nonamer_score = 0
                    if upstream_nonamer[0] == 'g': upstream_nonamer_score += 1
                    if upstream_nonamer[1] == 'g': upstream_nonamer_score += 1
                    if upstream_nonamer[2] == 't': upstream_nonamer_score += 1
                    if upstream_nonamer[3] == 't': upstream_nonamer_score += 1
                    if upstream_nonamer[4] == 't': upstream_nonamer_score += 1
                    if upstream_nonamer[5] == 't': upstream_nonamer_score += 1
                    if upstream_nonamer[6] == 't': upstream_nonamer_score += 1
                    if upstream_nonamer[7] == 'g': upstream_nonamer_score += 1
                    if upstream_nonamer[8] == 't': upstream_nonamer_score += 1
    
                    # Setting nonamer_score lower than 6 leads to more false positives and breaks results
                    if (upstream_heptamer[4:7] == 'gtg' and upstream_heptamer_score >= 2 
                        and upstream_nonamer_score >= 5):
                        sta_nt -= i
                        found = True
                        break
    
                gene     = nts[sta_nt:end_nt]
                gene_len = end_nt - sta_nt
                
                # Look for exact matches between discovered gene and known alleles
                # Allow for the possibility that multiple alleles have same sequence
                # Default assumption is that gene is not in database
    
                allele = 'Not in J ref db'
                for jallele, jseq in jseq_dict.items():
                    if gene == jseq:
                        if allele == 'Not in J ref db':
                            allele = jallele + ' ' + jtype_dict[jallele]
                        else:
                            allele = allele + ', ' + jallele + ' ' + jtype_dict[jallele]
                
                # Flag J genes that start before last V gene
                if last_v_nt and sta_nt < last_v_nt:
                    continue
                    annot += ' / located in V gene region'
    
                if found:
                    jout.write(f">{allele} {gene_len} nts: {sta_nt} - {end_nt}\n")
                    jout.write(gene)
                    jout.write("\n\n")
                else:
                    pseudogene_list.append(gene)
                    
    if pseudogenes:
        with open( pseudo_file, mode ) as pickle_file:
            pickle.dump( pseudogene_list, pickle_file ) 
        pickle_file.close()
        
    if last_v_nt:
        return last_v_nt
    else:
        print("last_v_nt =", last_v_nt)
                        
    jout.close()