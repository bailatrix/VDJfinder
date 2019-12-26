# VDJfinder
> A search tool for identifying V, D, and J gene segments in Ig and TCR loci. Optionally, supports identification of pseudogenes and scan of NCBI database.


#### TODO:
#### implement screening methods for necessary search parameters:
####    - locus_file
####    - locus_type
####    - gene_type
####    - pref_name
####    - custom_rules
#### implement support for outputing as .fasta file


## Example Usage
_This program has been developed for use in several platforms including command line execution, Jupyter Notebook, and our website. Below are examples of what running the search tool should look like on these platforms._

#### Command Line
    Method: python search( locus_file, locus_type, gene='ALL', custom_rules=False )

    Example: `python3 search( './data/input/IGH_locus.fasta', 'IGH')`

#### Jupyter Notebook
    Method: search( locus_file, locus_type, gene='ALL', custom_rules=False )

    Example: `python3 search( './data/input/IGH_locus.fasta', 'IGH')`

#### Original method
    Example:
    `python3 vdjfinder.py -f Loci/IGH_locus.fasta -v IMGT_human/IMGT/IGHV.fasta -d IMGT_human/IMGT/IGHD.fasta -j IMGT_human/IMGT/IGHJ.fasta`

#### Online
There is a web application in development as of January 2020. Once live, we will update with sample instructions.



## **Program Details**
> There is only one method intended to be public to the user, [search()](https://github.com/bailatrix/VDJfinder/blob/master/src/modules/search.py). The remainder of the program performs the search using a rigorously tested set of search criteria and prepares the result for the user. 


### Modules
search
* _Methods_
    * `search( locus_file, locus_type, gene='ALL', custom_rules=False )`

prepIO
* _Global values_
    * `all_ref_dbs`
    * NOTE: Nonamer matches >=5, heptamer matches >= 5, and no restriction on sum finds all J genes with no false positives

* _Methods_
    * `prep_output( gene_file, pseudogenes=False, pref_name=False, force=False )`
       - Private to v/d/j_gene_search() methods
       - preps output file to store search results
       - __Parameter(s)__:
             - `gene_file`: file that is being searched, naming convention uses this file name as base
             - `pseudogenes`: Boolean value, indicates whether to alter naming convention to include 'pseudogenes'
             - `pref_name`: allows out_name to be passed as argument
             - `force`: Boolean value, indicates whether to force overwrite in the event of duplicate file name
       - __Return__: 
             - file location+name and mode to use
         
    * `prep_database( locus_type, gene_type )`
       - Private to v/d/j_gene_search() methods
       - prepares local reference database as dictionaries
       - __Parameter(s)__:
             - `locus_type`: 'IGH', 'IGL', 'IGK', 'TRA', or 'TRB'
             - `gene_type`: 'V', 'D', or 'J'
       - __Return__:
             - dict of gene sequences and dict of gene types
    
    * `prep_frame( nt )`
       - Private to v_gene_search() method
       - Description TBD 
       - __Parameter(s)__:
             - `nt`: nucleotide sequence to operate on
       - __Return__:
             - amino acid frame

### Data
There is certain data which is vital for the search method's ability to remain accurate and is based on public data found on [The National Center for Biotechnology Information](https://www.ncbi.nlm.nih.gov/gene?Db=gene&Cmd=DetailsSearch&Term=3492) website. While we use a local version of this data in the program, we intend for this to be kept as up-to-date as possible. Please report this issue if you find this to not be the case. 

### Dependencies
__Python__ _(version 3.6)_
- Use: This program was written entirely in Python and cannot be executed outside of either an environment that supports Python or a third party environment designed to bridge this program with another platform (ie. web application). 
__Bio__ _(version 1.66)_
- Use: The Bio package is vital to this program's ability to read and write .fasta files, as well as the initial handling of each DNA sequence.
__Pickle__ _(version )_
- Use: During the search for V, D, and J genes, the program attempts to weed out false positives (AKA, "pseudogenes"). In order to continuously improve the accuracy of the program and further research this topic, when the program sorts out these genes, it collects the pseudogene(s) and pickles the collection for later investigation. While the goal of this collection is to increase knowledge of immunoglobulin genes across the board, this feature can be disabled by adding 'pseudogenes=False' as an argument to the initial call to search().


## **Contact**
> If you encounter any problems while using this program, please [report the bug](https://github.com/bailatrix/VDJfinder/issues) to the developers. Additionally, [contact](https://www.eloquenceintech.com/contact) the developers with any questions, comments, or problems running this tool.

This tool is a collaborative effort from researchers at the San Diego Supercomputer Center (SDSC) and Vanguard University.
- Bob Sinkovits, Ph.D. _Director of Scientific Computing Applications, SDSC_
- Bailey Passmore, Undergraduate Student, _Computational and Data Science Researcher, SDSC_
- _Additional names to be added_ 

![last updated](https://img.shields.io/github/last-commit/bailatrix/VDJfinder)
![repo size](https://img.shields.io/github/repo-size/bailatrix/VDJfinder)
[![issues](https://img.shields.io/github/issues/bailatrix/VDJfinder)](https://github.com/bailatrix/VDJfinder/issues)
![forks](https://img.shields.io/github/forks/bailatrix/VDJfinder?style=social)
![stars](https://img.shields.io/github/stars/bailatrix/VDJfinder?style=social)
[![Python 3.6](https://img.shields.io/badge/python-3.6-blue.svg)](https://www.python.org/downloads/release/python-360/)

---

Copyright 2019 San Diego Supercomputer Center