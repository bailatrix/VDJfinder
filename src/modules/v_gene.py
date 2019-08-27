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

# Search for v genes in given vgene_file based on given locus file
def v_search( locus,  vgene_file, last_v_nt=False, pseudogenes=False, overwrite=False ):
    pseudogene_list = []
    
    # Prepare output file and build local reference database
    v_file, mode = prep_IO.prep_output( vgene_file, force=overwrite )
    vout = open( v_file, mode )
    vseq_dict, vtype_dict = prep_IO.prep_database( vgene_file )
    
    if pseudogenes:
        pseudo_file, mode = prep_IO.prep_pseudo_file( vgene_file, force=overwrite )
    
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
    
        # -------------------- SEARCH FOR V GENES --------------------
    
        # Find matches for V-genes (1st conserved Cys to 2nd conserved Cys)
        # Max number of residues between Cys1-Trp       = 17
        # Max number of residues between Tpr-[IVLFCMA]  = 47
        # Max number of residues between [IVLFCMA]-Cys2 = 14
        # Allow for missing residues - 9 in Cys1-Trp, 9 in Trp-[IVLFCMA], 2 in [IVLFCMA]-Cys2
    
    
        for frame,offset in zip(aa_frame,[0,1,2]):
            for mcons in finditer('C[A-Z]{8,17}W[A-Z]{38,47}[IVLFCMA][A-Z]{12,14}C', frame):
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
    
                    # heptamer consensus = cacagtg
                    heptamer_score = 0
                    if heptamer[3] == 'a': heptamer_score += 1
                    if heptamer[4] == 'g': heptamer_score += 1
                    if heptamer[5] == 't': heptamer_score += 1
                    if heptamer[6] == 'g': heptamer_score += 1
                    
                    # nonamer consensus = acaaaaacc
                    nonamer_score = 0
                    if nonamer[0] == 'a': nonamer_score += 1
                    if nonamer[1] == 'c': nonamer_score += 1
                    if nonamer[2] == 'a': nonamer_score += 1
                    if nonamer[3] == 'a': nonamer_score += 1
                    if nonamer[4] == 'a': nonamer_score += 1
                    if nonamer[5] == 'a': nonamer_score += 1
                    if nonamer[6] == 'a': nonamer_score += 1
                    if nonamer[7] == 'c': nonamer_score += 1
                    if nonamer[8] == 'c': nonamer_score += 1
    
                    # Setting nonamer_score lower than 6 leads to more false positives and breaks results
                    if heptamer[0:3] == 'cac' and heptamer_score >= 2 and nonamer_score >= 6:
                        end_nt += i
                        post_cys2 = i
                        found = True
                        break
    
                # Set start of gene to 63 bases before first conserved Cys in IGH
                # (Use 66 bases for IGK and IGL)
                sta_nt  -= 63
                gene     = nts[sta_nt:end_nt]
                aa       = frame[sta_aa-21:end_aa+int(post_cys2/3)]
                gene_len = end_nt - sta_nt
                
                # Test V gene for stop codons
                if '*' in aa:
                    stop_codon = '(Contains stop codon)'
                else:
                    stop_codon = ''
    
                # Look for exact matches between discovered gene and known alleles
                # Allow for the possibility that multiple alleles have same sequence
                # Default assumption is that gene is not in database
                allele = 'Not in V ref db'
                for vallele, vseq in vseq_dict.items():
                    if gene == vseq:
                        if allele == 'Not in V ref db':
                            allele = vallele + ' ' + vtype_dict[vallele]
                        else:
                            allele = allele + ', ' + vallele + ' ' + vtype_dict[vallele]
                
                if found:                
                    vout.write(f">{allele} {gene_len} nts: {sta_nt} - {end_nt} {stop_codon}\n")
                    vout.write(gene)
                    vout.write("\n\n")
                    
                    last_v_nt = end_nt # Keep track of where last V gene found
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
    
    vout.close()