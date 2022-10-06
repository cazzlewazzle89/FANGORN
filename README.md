# FANGORN: Full-length Amplicons for the Next Generation Of rRNa analysis
### A quality-checked and publicly available database of full-length 16S-ITS-23S rRNA operon sequences

This repository makes available the scripts used to build the FANGORN databases described in this preprint (REMINDER TO LINK) and available for download [here](https://melbourne.figshare.com/articles/dataset/Fangorn_rrn_Database/20086916).  
FANGORN is envisaged as a tool to aid standardisation of 16S-ITS-23S rRNA analysis and allow comparison of results and, as such, building your own version would defeat the purpose.  

Please get in touch if you have any comments, issues, or suggestions for improvements.

I plan to update the database in line with each major RefSeq and GTDB release.

Note: If you want to build your own version using the NCBI taxonomy system, make sure you have the most up-to-date version of the taxonomy database. I do this using the commands described in the [TaxonKit manual](https://bioinf.shenwei.me/taxonkit/usage/#before-use).  
This does not download automatically as I often encounter an EOF error while extracting the `taxdump` tarball so it is less hassle to do it manually using the accomanying script `get_taxonkit_DB.sh` before running `refseq_farngorn.sh`.  
You do not need to do this when using the GTDB taxonomy system as the information is available on [the website](https://gtdb.ecogenomic.org/downloads) and will be automatically downloaded when running `gtdb_fangorn.sh`.

## Dependencies
Make sure these are in your $PATH

| Software  | Version Tested |
| --- | --- |
| [Barrnap](https://github.com/tseemann/barrnap) | 0.9 |
| [BBTools](https://jgi.doe.gov/data-and-tools/bbtools/) | 38.90  |
| [BEDTools](https://github.com/arq5x/bedtools2) | 2.30.0  |
| [csvtk](https://github.com/shenwei356/csvtk) | 0.24.0 |
| [R](https://www.r-project.org/) | 4.1.0  |
| [SeqKit](https://github.com/shenwei356/seqkit) | 2.2.0 |
| [TaxonKit](https://bioinf.shenwei.me/taxonkit/) | 0.8.0  |
| [VSEARCH](https://github.com/torognes/vsearch) | 2.17.1  |


### R Packages

| Package | Version Tested |
| --------|----------------|
| [Tidyverse](https://www.tidyverse.org/) | 1.3.1 |
| [Ape](https://cran.r-project.org/web/packages/ape/index.html) | 5.0 |

The conda environment that I used to build the database can be created using `conda env create -f fangorn.yml` (the yml is provded in this repository).  
This will create an environment called `fangorn` which can be loaded with `conda activate fangorn`.  
This environment contains everything except `csvtk` and `Seqkit`. These are available through bioconda but I used preexisting non-conda installations.

## Usage

`refseq_fangorn.sh` and `gtdb_fangorn.sh` are used to build databases using RefSeq and GTDB taxonomy respectively.  
By default they will use:  
- 50 parallel processes to download assemblies and annotations and edit seqid headers  
- 50 threads for VSEARCH clustering  
- 20 threads for identifying rRNA genes in genomes without accompanying annotation using Barrnap  

These settings, along with VSEARCH clustering identity (default: 99.9%) and output directory (default: `./Fangorn_GTDB/` or `./Fangorn_RefSeq/`) can be changed by editing the variables at the beginning of the relevant script.  

Both scripts mentioned above will do the majority of the work but will call R scripts to identify potential _rrn_ operons, perform quality checking (the details of which can be found in the FANGRON publication), and join _rrn_ operon identifiers to the source genome taxonomy.
R scripts are also used by the GTDB database builder to extract assembly information and identify genomes which are missing annotation.  

Also included in the `Scripts/` folder are scripts to:  
* collate the files necessary to generate the manuscripts figures and statistics  
* calculate average nucleotide identity (ANI) between all pairwise combinations of complete genomes used to construct both GTDB and RefSeq databases  
* train a Qiime2 Naive Bayes feature classifier to assign taxonomy to amplicons (generated using all tested primer pairs) using each database - I will make these available shortly

## To-do  

:muscle: Train NB classifier for each primer set based on each database and taxonomy system  
:microscope: Compare phylogenetic resolution of V4, 16S and _rrn_  (is full 16S or _rrn_ more representative of whole genome identity that V4?)  
:palm_tree: Build phylogenetic tree for SEPP  