# FANGORN: Full-length Amplicons for the Next Generation Of rRNa analysis
### A quality-checked and publicly available database of full-length 16S-ITS-23S rRNA operon sequences

This repository makes available the scripts used to build the FANGORN databases described in this preprint (REMINDER TO LINK) and available for download [here](https://melbourne.figshare.com/account/projects/119793/articles/19530565).
FANGORN was envisaged as a tool to aid standardisation of 16S-ITS-23S rRNA analysis and allow comparison of results and, as such, building your own version would defeat the purpose.
Please get in touch if you have any comments, issues, or suggestions for improvements.

I plan to update the database with each new GTDB release.

Note: If you want to build your own version using the NCBI taxonomy system, make sure you have the most up-to-date version of the taxonomy database. I do this using the commands described in the [TaxonKit manual](https://bioinf.shenwei.me/taxonkit/usage/#before-use) described below.
```bash
wget -c ftp://ftp.ncbi.nih.gov/pub/taxonomy/taxdump.tar.gz
tar -zxvf taxdump.tar.gz

mkdir -p TaxonKit/
mv names.dmp nodes.dmp delnodes.dmp merged.dmp TaxonKit/
rm -f citations.dmp division.dmp gc.prt gencode.dmp readme.txt taxdump.tar.gz
```
You do not need to do this when using the GTDB taxonomy system as the information is available on [the website](https://gtdb.ecogenomic.org/downloads) and will be automatically downloaded in the scripts provided.

## Dependencies
Make sure these are in your $PATH

| Software  | Version Tested |
| --- | --- |
| [Barrnap](https://github.com/tseemann/barrnap) | 0.9 |
| [BBTools](https://jgi.doe.gov/data-and-tools/bbtools/) | 38.90  |
| [BEDTools](https://github.com/arq5x/bedtools2) | 2.30.0  |
| [R](https://www.r-project.org/) | 4.1.0  |
| [TaxonKit](https://bioinf.shenwei.me/taxonkit/) | 0.8.0  |
| [VSEARCH](https://github.com/torognes/vsearch) | 2.17.1  |


### R Packages

| Package | Version Tested |
| --------|----------------|
| [Tidyverse](https://www.tidyverse.org/) | 1.3.1 |
| [Ape](https://cran.r-project.org/web/packages/ape/index.html) | 5.0 |

## Usage

`refseq_fangorn.sh` and `gtdb_fangorn.sh` are used to build databases using RefSeq and GTDB taxonomy respectively.
By default they will use 50 parallel processes to download assemblies and annotations and edit seqid headers, 50 threads for VSEARCH clustering, and 20 threads for identifying rRNA genes in genomes without accompanying annotation using Barrnap. These settings, along with Vsearch clustering identity (default: 99.9%) and output directory (default: ${PWD}/Fangorn_[GTDB/RefSeq]) can be changed by editing the relevant values in lines 2-7 of the script. 