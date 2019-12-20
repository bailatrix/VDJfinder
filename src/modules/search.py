from modules import v_gene
from modules import d_gene
from modules import j_gene

############################
# TODO:
#    - implement support for searching multiple locus_files at a time?
############################


# search_db():
#    - Custom search through db for pseudogenes
# Parameters:
#    - 
# Return:
#    - 
def search_db( ):
    
    return 0

# This method is intended to be THE public search called by users
# Notes: 
#    - Search logic is specific to each locus+gene type pair
#         - ie. given locus = 'IGH' and gene = 'V', search will follow prescribed 'IGHV' rules
#    - D and J gene locations must be checked against the location of the last known V nucleotide, last_v_nt
#        - If found before the location of last_v_nt, these genes are determined to be erroneous and misidentified
# search():
#    - Takes user call and initiates corresponding searches
# Parameters:
#    - locus_file: file to search
#    - locus_type: 'IGH', 'IGL', 'IGK', 'TRA', or 'TRB'
#    - gene: Defaults to 'ALL', otherwise can be 'V', 'D', or 'J'
#    - custom_rules: Boolean value, defaults to False. If True, indicates user wants to search file with custom search values
# Return:
#    - No return value, but search programs called create respective '_found' files for genes
def search( locus_file, locus_type, gene='ALL', custom_rules=False ):
    from sys import exc_info
    
    try:
        
        locus_type = locus_type.upper()
        gene = gene.upper()
        
        # Custom rules must be handled by separate search functions to handle IO properly
        if custom_rules:
            if (gene == 'ALL') and (locus_type in ['IGH', 'TRB']):
                
                last_v_nt = v_gene.custom_v_search( locus_file, locus_type )
                d_gene.custom_d_search( locus_file, locus_type, last_v_nt )
                j_gene.custom_j_search( locus_file, locus_type, last_v_nt )
        
            elif (gene == 'ALL') and (locus_type in ['IGK', 'IGL', 'TRA']):
        
                last_v_nt = v_gene.custom_v_search( locus_file, locus_type )
                j_gene.custom_j_search( locus_file, locus_type, last_v_nt )        
        
        else:
            # Given standard search rules
            if (gene == 'ALL') and (locus_type in ['IGH', 'TRB']):
                
                last_v_nt = v_gene.v_search( locus_file, locus_type )
                d_gene.d_search( locus_file, locus_type, last_v_nt )
                j_gene.j_search( locus_file, locus_type, last_v_nt )
                
            elif (gene == 'ALL') and (locus_type in ['IGK', 'IGL', 'TRA']):
                
                last_v_nt = v_gene.v_search( locus_file, locus_type )
                j_gene.j_search( locus_file, locus_type, last_v_nt )        
    
    # handle incorrect number of arguments passed
    except TypeError:
        print("Invalid input. search() takes at least 2 string inputs as arguments: locus_file and locus_type.")
        raise TypeError
    
    # handle invalid arguments passed
    except ValueError:
        print("Invalid input. search() takes at least 2 string inputs as arguments: locus_file and locus_type.")
        raise ValueError
    
    # handle unknown error
    except:
        print("Unexpected error:", exc_info()[0])
        raise
        
if __name__ == '__main__':
    search(sys.argv)