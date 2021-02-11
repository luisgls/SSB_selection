## SSB-dN/dS Method to detect positive and negative selection in cancer
### *Luis Zapata*

SSB dN/dS is a tool that calculates the ratio of nonsynonymous to synonymous mutations (dN/dS) for genes using annotated variant data. Given that some mutational processes are more common than others, SSB-dN/dS uses a context correction method that normalizes for this bias (Somatic Substitution Bias, SSB). The tool and the results are described in  the manuscript:

["Negative Selection on tumour evolution acts on essential cellular functions and the immunopeptidome"](https://doi.org/10.1186/s13059-018-1434-0).

## Installation

### Dependencies
bedtools 2.26.0 
https://github.com/arq5x/bedtools2

R-3.3.3 or higher

Python (to install the synapse client)

R library tidyr

perl 5

variant effect predictor v89 or higher

GNU command line tools



#### Important Notes
- earlier versions of bedtools will not work
- tab encoding should be \t (might be a problem for windows/OSX versions)
- genome file is a two column file specifying the fasta id and the length of the sequence (see how to obtain it at the bottom)
- Restrict your input dataset to chromosomes 1-22 and X and Y. Remove the rest.

#### To install first clone the tool

```{bash}
git clone https://github.com/luisgls/SSB_selection.git
```
#### Install synapse client  *need pip from python
```{bash}
pip install synapseclient
```

#### Download zipped data files from synapse syn11681952
```{bash}
synapse get -r syn11681952
```

#### Go to the tool directory and
  1.1) Create data folder within the cloned folder
  
  1.2) Unzip data files into data folder (Data.zip for HG19, and Data2.zip for GRCh38)
  
  1.3) Unzip Example_files into example directory 
  
  1.4) Open run_negDriver script and specify the location of the genome file (e.g. hg19.genome) and the fasta file (e.g. hg19.fasta) 
  
  1.5) You need to also mody the BASEDIR parameter where your currrent instalation of the scripts are.
  
  1.5) chmod +x run_negDriver
  
2) run ./run_negDriver to see help

## Input file
The input file is the standard output of variant effect predictor using the following command line (by providing to vep the ensembl default input file format)

```{bash}
perl variant_effect_predictor.pl -i input -o input.annotated --cache --all_refseq --assembly GRCh37 --pick --symbol --no_stats --fasta genome.fasta
```

If you want to filter putative germline variants use the option --plugin ExAC when running VEP.

Example input files can be found on synapse: ID syn11681983

#### Important points before running
a) No header needed for input VEP file

b) VEP annotated first column must be in the format (chr_pos_ref/alt)

c) VEP annotated file must only have chromosomes that are 1,2,3,4...22 or uppercase X,Y

d) After you run vep with the option for ExAC frequencies, it would be necessary to remove all variants present in more than 0.1 percent of the population. You could apply the filter using:
```{bash}
filter_vep -i input.annotated -f "ExAC_AF < 0.1 or not ExAC_AF" --ontology --filter "Consequence is coding_sequence_variant" 
```

e) Be sure that you are using the GNU command line if you are running in a MacOS (https://www.topbug.net/blog/2013/04/14/install-and-use-gnu-command-line-tools-in-mac-os-x/)

f) Add dependencies to your path for easy running or hardcode the scripts

g) The genome file used in SSB is a two column file that contains the info of the name of the fasta id (column 1) and the length of that sequence (column 2).

h) The UNIX system used should be able to recognize \t as a tab separator

NOTES:
This version is a simplified and easy-to-use version of the tool developed in the manuscript "Negative Selection on tumour evolution acts on essential cellular functions and the immunopeptidome" (https://doi.org/10.1186/s13059-018-1434-0). For the reproducibility of the results presented in the manuscript, please contact the authors directly.

The tool to run the analysis of the immunopeptidome is called [SOPRANO](https://github.com/luisgls/SOPRANO)

#### Run the tool for other reference assemblies
To run the tool with version 38 of the human genome, simply update the path of you GENOME and FASTA file in the run_negdriver script.
Also, use the Data file provided in synapse labeled as Data2.zip. It contains updated transcript information.

## Genomes
To get hg19 fasta genome, you can download it from UCSC:

```{bash
wget http://hgdownload.cse.ucsc.edu/goldenPath/hg19/bigZips/hg19.fa.gz

wget http://hgdownload.cse.ucsc.edu/goldenPath/hg19/bigZips/hg19.chrom.sizes

```

Remember to modify the fasta file to contain only uppercase letters.
