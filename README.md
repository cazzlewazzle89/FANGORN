# FANGORN: Full-length Amplicons for the Next Generation Of rRNa analysis
### A quality-checked and publicly available database of full-length 16S-ITS-23S rRNA operon sequences

This repository makes available the scripts used to build the FANGORN databases described this this preprint and available for download here (REMINDER TO LINK).  
FANGORN was envisaged as a tool to aid standardisation of 16S-ITS-23S rRNA analysis and allow comparison of results and, as such, building your own version would defeat the purpose.  
Please get in touch if you have any comments, issues, or suggestions for improvements.

I plan to update the database with each new GTDB release. For the RefSeq/NCBI version - every 12 months seems logical.

Note: If you want to build your own version using the NCBI taxonomy system, make sure you have the most up-to-date version of the taxonomy database. I do this using the commands described in the [TaxonKit manual](https://bioinf.shenwei.me/taxonkit/usage/#before-use).

```
wget -c ftp://ftp.ncbi.nih.gov/pub/taxonomy/taxdump.tar.gz 
tar -zxvf taxdump.tar.gz

mkdir -p TaxonKit/
mv names.dmp nodes.dmp delnodes.dmp merged.dmp TaxonKit/
rm -f citations.dmp division.dmp gc.prt gencode.dmp readme.txt taxdump.tar.gz
```

## Dependencies  
Make sure these are in your $PATH  

| Software  | Version Tested |
| --- | --- |
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
