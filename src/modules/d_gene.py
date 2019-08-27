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

# Search for d genes in given dgene_file based on given locus file
def d_search( locus, dgene_file, last_v_nt=False, pseudogenes=False, overwrite=False ):    
    pseudogene_list = []
    
    # Prepare output file and build local reference database
    d_file, mode = prep_IO.prep_output( dgene_file, force=overwrite )
    dout = open( d_file, mode )
    dseq_dict, dtype_dict = prep_IO.prep_database( dgene_file )
    
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
        
        #for dseq in re.finditer('cac[acgt]gtg[acgt]{10,32}cac[acgt]gtg', nts):
        #w/ hept>=2, non>=4 16 true pos, 0 false pos
    
        #for dseq in re.finditer('[acgt]ac[acgt]gtg[acgt]{10,32}cac[acgt]gtg', nts):
        # w/ hept>=2, non>=4 22 true pos, 2 false pos
    
        for dseq in finditer('[acgt]ac[acgt]gtg[acgt]{10,37}cac[acgt]gtg', nts):
            sta_nt = dseq.span()[0]
            end_nt = dseq.span()[1]
            gene = nts[sta_nt+7:end_nt-7]
            spacer = nts[end_nt:end_nt+12]
            gene_len = end_nt - sta_nt - 14
            
            # upstream heptamer consensus = cactgtg
            upstream_heptamer = nts[sta_nt:sta_nt+7]
            upstream_heptamer_score = 0
            if upstream_heptamer[0] == 'c': upstream_heptamer_score += 1
            if upstream_heptamer[1] == 'a': upstream_heptamer_score += 1
            if upstream_heptamer[2] == 'c': upstream_heptamer_score += 1
            if upstream_heptamer[3] == 't': upstream_heptamer_score += 1
    
            # upstream nonamer consensus = tgtttttgt
            upstream_nonamer = nts[sta_nt-12-9:sta_nt-12]
            upstream_nonamer_score = 0
            if upstream_nonamer[0] == 't': upstream_nonamer_score += 1
            if upstream_nonamer[1] == 'g': upstream_nonamer_score += 1
            if upstream_nonamer[2] == 't': upstream_nonamer_score += 1
            if upstream_nonamer[3] == 't': upstream_nonamer_score += 1
            if upstream_nonamer[4] == 't': upstream_nonamer_score += 1
            if upstream_nonamer[5] == 't': upstream_nonamer_score += 1
            if upstream_nonamer[6] == 't': upstream_nonamer_score += 1
            if upstream_nonamer[7] == 'g': upstream_nonamer_score += 1
            if upstream_nonamer[8] == 'g': upstream_nonamer_score += 1 # RSS ERROR
    
            # downstream heptamer consensus = cacagtg
            downstream_heptamer = nts[end_nt-7:end_nt]
            downstream_heptamer_score = 0
            if downstream_heptamer[3] == 'a': downstream_heptamer_score += 1
            if downstream_heptamer[4] == 'g': downstream_heptamer_score += 1
            if downstream_heptamer[5] == 't': downstream_heptamer_score += 1
            if downstream_heptamer[6] == 'g': downstream_heptamer_score += 1
    
            # downstream nonamer consensus = acaaaaacc
            downstream_nonamer = nts[end_nt+12:end_nt+12+9]
            downstream_nonamer_score = 0
            if downstream_nonamer[0] == 'a': downstream_nonamer_score += 1
            if downstream_nonamer[1] == 'c': downstream_nonamer_score += 1
            if downstream_nonamer[2] == 'a': downstream_nonamer_score += 1
            if downstream_nonamer[3] == 'a': downstream_nonamer_score += 1
            if downstream_nonamer[4] == 'a': downstream_nonamer_score += 1
            if downstream_nonamer[5] == 'a': downstream_nonamer_score += 1
            if downstream_nonamer[6] == 'a': downstream_nonamer_score += 1
            if downstream_nonamer[7] == 'c': downstream_nonamer_score += 1
            if downstream_nonamer[8] == 'c': downstream_nonamer_score += 1
    
    
            # Look for exact matches between discovered gene and known alleles
            # Allow for the possibility that multiple alleles have same sequence
            # Default assumption is that gene is not in database
    
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
                allele += ' / located in V gene region'
    
            if (upstream_heptamer_score >= 2 and upstream_nonamer_score >= 4 and 
                upstream_heptamer_score >= 2 and upstream_nonamer_score >= 4):
                
                dout.write(f">{allele} {gene_len} nts: {sta_nt} - {end_nt}\n")
                dout.write(gene)
                dout.write("\n\n")  
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
                
    dout.close()