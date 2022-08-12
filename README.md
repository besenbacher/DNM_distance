# DNM_distance

## Requirements
[py2bit](https://github.com/deeptools/py2bit) needs to be installed for the python scripts to work. This can be done using pip:
```
pip install git+https://github.com/dpryan79/py2bit
```
For the R plotting script to work the `tidyverse` and `Hmisc` packages needs to be installed.


## Input files
The de novo mutations should be in a txt file with the following columns: chrom, pos, ref, alt, sample_id. Where sample_id is an identifier of which proband the mutation occured in.
The program also requires a bed-file listing which regions are sufficiently covered by the sequencing. Both files needs to be sorted first by chromosome (alphabetically) and then position (numerically). 

If no bed-file is available a bed-file listing all chromosomes as covered can be created from a list of chromosome lengths like this:
```
cat chromSizes_hg38.txt | awk -v OFS="\t" '{print $1,"0",$2}' | sort -k1,1 -k2,2n > all_chromosomes_hg38.bed
```
The program also requires the a 2bit files with the reference sequence. Such a file can be downloaded from:
`https://hgdownload.cse.ucsc.edu/goldenpath/{genome}/bigZips/{genome}.2bit` where `{genome}` is a valid UCSC genome assembly name (fx. "hg38").

## Example
The DNM_distance script can be run on the example input data from the test_data directory like this:
```
python src/DNM_distance.py test_data/decode2019_mutations_chr21_chr22_hg38.gor hg38.2bit --bed test_data/all_chromosomes_chr21_chr22_hg38.bed -n100 > test_output_chr21_chr22_n100.txt
```
Make plot:
```
./R_scripts/plot_DNM_distances.R test_chr21_chr22_n100.txt test_chr21_chr22_n100_DNM_distances.pdf
```
To analyse all pairwise mutations:
```
python src/DNM_distance/all_pair_distances.py test_data/decode2019_mutations_chr21_chr22_hg38.gor| gzip > all_pair.txt.gz
```
Make plot:
```
./R_scripts/plot_all_pair_results.R all_pair.txt.gz 2976 all_pair_results.pdf
```
The number is 2976 is the number of probands in the trio data and should be replaced with the right number when the script is run on another dataset.
