def search_file( ):
    
    return 0

def search_db( ):
    
    return 0

def search( ):
    
    # Read fasta sequence (ntf = forward , ntr = reverse)
    for nt in SeqIO.parse(locus_file, "fasta"):
    
        if read_dir == 'reverse':
            nt = nt.reverse_complement()
    
        nt_len = len(nt.seq)
        nts = str(nt.seq).lower()
    
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
    
        if vgene_file:
            last_v_nt = v_gene_search(vgene_out, aa_frame, last_v_nt, vseq_dict, v_motif, v_gap, 
                                      v_heptamer_consensus, v_nonamer_consensus, v_nt_before_cys1, 
                                      v_min_heptamer_match, v_min_nonamer_match, v_min_total_match)
    
        if dgene_file:
            d_gene_search(dgene_out, nts, last_v_nt, dseq_dict, d_motif, 
                          d_upstream_heptamer_consensus, d_upstream_nonamer_consensus, 
                          d_downstream_heptamer_consensus, d_downstream_nonamer_consensus, 
                          d_min_upstream_heptamer_match, d_min_upstream_nonamer_match, 
                          d_min_downstream_heptamer_match, d_min_downstream_nonamer_match,
                          d_min_total_match)
    
        if jgene_file:
            j_gene_search(jgene_out, aa_frame, last_v_nt, jseq_dict,
                          j_motif, j_gap, j_heptamer_consensus, j_nonamer_consensus, 
                          j_min_heptamer_match, j_min_nonamer_match, j_min_total_match)