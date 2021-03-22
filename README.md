# VCF Tempus Parser
This script annotates a VCF file and outputs a CSV or TSV

## Installation
Clone the Tempus repo directory from GitHub:

```
ticks
git clone https://github.com/jamesthesnake/Tempus_repo.git
cd Tempus_repo
```
Install requirements.txt:

```
pip3 install -r requirements.txt
You should have python 3.7+ installed
```

Install python3 and install the requests package
pip install requirements.txt
this is the use of the program

## Usage
Run:
```
var_tempus_annotate.py [-h] vcf outputfile

positional arguments:
  vcf         VCF file location, put the input file here
  outputfile  output file name, put the name of the output file you would like

optional arguments:
  -h, --help  show this help message and exit
```
## Output glossary
|Key | Description|
|----|------------|
|Chr | Chromosome variant is locationed on|
|Position | Variants' location on the chromosome|
|Ref | Reference Allele|
|Alt | Alternate Allele|
|DP | Total depth|
|Annotation | Type of variant e.g. snp, insertion, del|
|SequenceDepth | Total read depth at variant location|
|VarsReads | Total number of reads representing the variant i.e. Variant Supporting Reads|
|PerReads | Percent of reads representing the variant|
|Consequences | What the consquence of the variant is |
|AlleFrequency | Allele Frequency according to ExAC Browser - Harvard|
|Genes | Affected gene according to ExAC Browser - Harvard|
######
