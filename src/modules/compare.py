#########################################
# Setting up comparison logic for genes #
#########################################


# Intended to be private to the compare_dict() function
# - compare():
#     - returns indicator based on if the first argument is longer or shorter than the second
#     - designed to support Python 'str1 in str2' search method

def compare( gene, seq ):
    if len(gene) > len(seq):
        return 1
    elif len(gene) < len(seq):
        return -1
    return 0


# Intended to be private to the compare_dict() function
# - compare_substring():
#     - returns dictionary of partial matches found for given gene,seq argument
#     - designed to continue search for a match in seq_dict by looking up all possible substrings of gene

def compare_substring( gene, seq, allele ):
    partials = {}
    
    for i in range(len(gene)):
        for j in range(1, len(gene)+1):
            partial = gene[i:-j:1]
            
            if partial in seq:
                match = round( len(partial)/len(seq), 3 )
                
                if (allele not in partials) and (match > 0.1):
                    partials[allele] = match
                elif allele in partials and partials[allele] < match:
                    partials[allele] = match
                    
                break
    return partials


# Intended to be private to the compare_dict() function
# - compare_partials():
#     - returns dictionary of all partial matches for gene across seq_dict
#     - designed to behave like compare_dict() but perform more intricate search for match

def compare_partials( gene, seq_dict ):
    matches = {}
    
    for allele, seq in seq_dict.items():
        partials = compare_substring(gene, seq, allele)
        
        if partials:
            matches.update(partials)
            
    return matches


# Called by the main program
# - compare_dict():
#     - returns dictionary of matches for given gene in seq_dict
#     - keys are the allele argument of seq, value is match percent

def compare_dict( gene, seq_dict ):
    matches = {}
    
    for allele, seq in seq_dict.items():
        comp = compare(gene, seq)
        match_perc = round(len(gene)/len(seq), 3)
        
        if comp > 0:
            if seq in gene:
                matches[allele] = match_perc
        elif comp < 0:
            if gene in seq:
                matches[allele] = match_perc
        else:
            if gene == seq:
                matches[allele] = match_perc
                
    if matches:
        return matches
    else:
        return compare_partials( gene, seq_dict )

# Called by the main program
# - pseudo_matches():

def pseudo_matches( pseudo_list, seq_dict ):
    matches = {}

    for gene in pseudo_list:
        matches[gene] = compare_dict( gene, seq_dict )
    
    return matches
    
# Called by the main program
# - top_matches():

def top_matches( partials_dict, minimum=0.25 ):
    top_matches = []

    for gene,matches in partials_dict.items():
        for allele,perc in matches.items():
            if perc >= minimum:
                top_matches.append( (allele,perc,gene) )
    return top_matches