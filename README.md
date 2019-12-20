# VDJfinder

A search tool for identifying V, D and J gene segments in Ig and TCR loci. Optionally, includes identification of pseudogenes.

### Example Usage
-----------------

#### Command Line
Method: python search( locus_file, locus_type, gene='ALL', custom_rules=False )

Example: python3 search( './data/input/IGH_locus.fasta', 'IGH')

#### Jupyter Notebook
Method: search( locus_file, locus_type, gene='ALL', custom_rules=False )

Example: python3 search( './data/input/IGH_locus.fasta', 'IGH')


#### Online
There is a web application in development as of January 2020. Once live, we will update with sample instructions.



Example: python3 vdjfinder.py -f Loci/IGH_locus.fasta -v IMGT_human/IMGT/IGHV.fasta -d IMGT_human/IMGT/IGHD.fasta -j IMGT_human/IMGT/IGHJ.fasta

### Program Details
---------------------
There is only one method intended to be public to the user, [search()](../src/modules/search.py). The remainder of the program performs the search using a rigorously tested set of search criteria and prepares the result for the user. 

#### Modules
search
⋅⋅* search( locus_file, locus_type, gene='ALL', custom_rules=False )

prepIO
.._Global values_
....* all_ref_dbs: 

.._Methods_
....* prep_output( gene_file, pseudogenes=False, pref_name=False, force=False ):  
....* prep_database( locus_type, gene_type ):
....* prep_frame( nt ):

##### Data


#### Authors
This tool is a collaborative effort from researchers at the San Diego Supercomputer Center (SDSC) and Vanguard University.
Bob Sinkovits, Ph.D. _Director of Scientific Computing Applications, SDSC_
Bailey Passmore, Undergraduate Student, _Computational and Data Science Researcher, SDSC_
Additional

### Contact
If you encounter any problems while using this program, please [report the bug](https://github.com/bailatrix/VDJfinder/issues) to the developer! Contact us with any questions, comments, or problems running this tool.

Updated: [![Updated](https://img.shields.io/github/last-commit/bailatrix/VDJfinder)]

[![repo size](https://img.shields.io/github/repo-size/bailatrix/VDJfinder)]
[![issues](https://img.shields.io/github/issues/bailatrix/VDJfinder)](https://github.com/bailatrix/VDJfinder/issues)
[![forks](https://img.shields.io/github/forks/bailatrix/VDJfinder?style=social)]
[![stars](https://img.shields.io/github/stars/bailatrix/VDJfinder?style=social)]

Copyright [2019] [San Diego Supercomputer Center]

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

    http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.