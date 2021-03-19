# VCF Parser

This script parses a VCF file and outputs a table in tab-separated values (TSV) form.  
The table contains the following columns:  

|Key | Description|
|----|------------|
|Chr | Chromosome variant is locationed on|
|Position | Variants' location on the chromosome|
|Ref | Reference Allele|
|Alt | Alternate Allele|
|Annotation | Type of variant e.g. snp, insertion, del|
|SeqDepth | Total read depth at variant location|
|VSReads | Total number of reads representing the variant i.e. Variant Supporting Reads|
|PerReads | Percent of reads representing the variant|
|Consequences | What the consquence of the variant is |

|AlleFrequency | Allele Frequency according to ExAC Browser - Harvard|
|Genes | Affected gene according to ExAC Browser - Harvard|
Clone the Tempus Repo directory from Github
######

Cd into Tempus_repo from your terminal shell

in the terminal console type: pip install -r requirements.txt
or
pip3 install -r requirements.txt
You should have python 3.7+ installed

Install python3 and install the requests package
pip install requirements.txt

This is the programi


usage: var_tempus_annotate.py [-h] vcf outputfile

positional arguments:
  vcf         VCF file location, put the input file here
  outputfile  output file name, put the name of the output file you would like

optional arguments:
  -h, --help  show this help message and exit



