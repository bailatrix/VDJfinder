# --- Default model parameters ---

# This is where we set the default parameters for the the V, D and J
# gene models. This includes motifs that are used to find the genes,
# heptamer and nonamer consensus sequences and required number of
# matches for nonamers, heptamers or combination. Note that we ALWAYS
# require the three heptamer bases adjacent to the V, D, or J gene to
# exactly match the consensus.

# -- IGHV --
# Motif description: Searching for pair of conserved Cys with intervening conserved Trp and [IVLFCMA]
# Motif description: 8-17 aa between Cys1-Trp, 38-47 aa between Trp-[IVLFCMA], 12-14 aa between [IVLFCMA]-Cys2
# Parameter notes: Nonamer matches < 6 or heptamer matches < 5 lead to more false positives
# Parameter notes: Nonamer matches > 6 or heptamer matches > 5 lose real genes
# Parameter notes: Total matches >= 13 filters ORFs without losing real genes
ighv_motif = 'C[A-Z]{8,17}W[A-Z]{38,47}[IVLFCMA][A-Z]{12,14}C'
ighv_heptamer_consensus = 'cacagtg'
ighv_nonamer_consensus = 'acaaaaacc'
ighv_nt_before_cys1 = 63
ighv_min_heptamer_match = 5
ighv_min_nonamer_match = 6
ighv_min_total_match = 13

# -- IGHD --
# Motif description: Due to short length and lack of highly conserved regions, motif needs to include heptamers
# Motif description: Look for 10-37 bases flanked by three conserved nt on each side plus additional conserved nt

ighd_motif = '[acgt]ac[acgt]gtg[actg]{10,37}cac[acgt]g[actg]{2}'
ighd_upstream_heptamer_consensus = 'cactgtg'
ighd_upstream_nonamer_consensus = 'tgtttttgg'   #RSS question of whether last base is g or t
ighd_downstream_heptamer_consensus = 'cacagtg'
ighd_downstream_nonamer_consensus = 'acaaaaacc'
ighd_min_upstream_heptamer_match = 5
ighd_min_upstream_nonamer_match = 4
ighd_min_downstream_heptamer_match = 5
ighd_min_downstream_nonamer_match = 4
ighd_min_total_match = 22


# -- IGHJ --
# Motif description: Trp and Ser-Ser separated by exactly 8 aa
# Parameter notes: Nonamer matches >=5, heptamer matches >= 5 and no restriction on sum finds all J genes with no false positives
ighj_motif = 'W[A-Z]{8}SS'
ighj_heptamer_consensus = 'cactgtg'
ighj_nonamer_consensus = 'ggtttttgt'
ighj_min_heptamer_match = 5
ighj_min_nonamer_match = 5
ighj_min_total_match = 0