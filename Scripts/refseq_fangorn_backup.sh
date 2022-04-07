#!/bin/sh
VAR_PARALLEL_DL=50 # number of parallel processes to run when downloading assemblies & annotations
VAR_PARALLEL_EDIT=50 # number of parallel processes to run when editing seqid headers
VAR_THREADS_VSEARCH=50 # number of threads to use then clustering operon sequences with vsearch
VAR_CLUSTERID_VSEARCH=0.999 # vsearch clustering identity
VAR_DB_TAXONDB="/home/cwwalsh/Databases/Taxonkit_DB/" # path to Taxonkit database
VAR_OUTPUT_DIRECTORY=${PWD}/Fangorn_RefSeq # specify output directory name for databases

# make output directory and move shell to run from there
mkdir ${VAR_OUTPUT_DIRECTORY}
cd ${VAR_OUTPUT_DIRECTORY}

# make directories to store genomes used for building the database
mkdir -p Genomes_Complete/
mkdir -p Genomes_Incomplete/

# get information for the lastest version of complete prokaryotic genome assemblies from RefSeq
wget https://ftp.ncbi.nlm.nih.gov/genomes/refseq/bacteria/assembly_summary.txt \
    -O RefSeq/bacteria_assembly_summary.txt
wget https://ftp.ncbi.nlm.nih.gov/genomes/refseq/archaea/assembly_summary.txt \
    -O RefSeq/archaea_assembly_summary.txt

cat RefSeq/archaea_assembly_summary.txt RefSeq/bacteria_assembly_summary.txt > RefSeq/prokaryote_assembly_summary.txt
rm -f RefSeq/archaea_assembly_summary.txt RefSeq/bacteria_assembly_summary.txt
awk -F '\t' '$11=="latest"' \
    RefSeq/prokaryote_assembly_summary.txt > RefSeq/prokaryote_assembly_summary_latest.txt
awk -F '\t' '$12=="Complete Genome"' \
    RefSeq/prokaryote_assembly_summary_latest.txt > RefSeq/prokaryote_assembly_summary_latestcomplete.txt
awk -F '\t' '$12!="Complete Genome"' \
    RefSeq/prokaryote_assembly_summary_latest.txt > RefSeq/prokaryote_assembly_summary_latestincomplete.txt

# extract ftp directory paths for each assembly
cut -f 20 RefSeq/prokaryote_assembly_summary_latestcomplete.txt > RefSeq/ftpdirpathscomplete
cut -f 20 RefSeq/prokaryote_assembly_summary_latestincomplete.txt > RefSeq/ftpdirpathsincomplete

# compress assembly information
# retaining to extract taxids and keep information about database version
gzip RefSeq/prokaryote_assembly_summary.txt 
gzip RefSeq/prokaryote_assembly_summary_latest.txt 
gzip RefSeq/prokaryote_assembly_summary_latestcomplete.txt 
gzip RefSeq/prokaryote_assembly_summary_latestincomplete.txt

# download nucleotide sequences
awk 'BEGIN{FS=OFS="/";filesuffix="genomic.fna.gz"}{ftpdir=$0;asm=$10;file=asm"_"filesuffix;print "wget -P RefSeq/Genomes_Complete/ "ftpdir,file}' \
    RefSeq/ftpdirpathscomplete > RefSeq/ftpdownload
parallel -a RefSeq/ftpdownload -j ${VAR_PARALLEL_DL}
rm -f RefSeq/ftpdownload

awk 'BEGIN{FS=OFS="/";filesuffix="genomic.fna.gz"}{ftpdir=$0;asm=$10;file=asm"_"filesuffix;print "wget -P RefSeq/Genomes_Incomplete/"ftpdir,file}' \
    RefSeq/ftpdirpathsincomplete > RefSeq/ftpdownload
parallel -a RefSeq/ftpdownload -j ${VAR_PARALLEL_DL}
rm -f RefSeq/ftpdownload

# simplify assembly filenames
find RefSeq/Genomes_Complete/ -type f -name *.fna.gz | awk -F '.[0-9]_' '{print "mv "$0" "$1".fna.gz"}' | bash
find RefSeq/Genomes_Incomplete/ -type f -name *.fna.gz | awk -F '.[0-9]_' '{print "mv "$0" "$1".fna.gz"}' | bash

# edit fasta headers to contain the RefSeq assembly accession
for i in $(find RefSeq/Genomes_Complete/ -type f -name *.fna.gz | sed 's/RefSeq\/Genomes_Complete\///' | sed 's/.fna.gz//')
do
    echo mv RefSeq/Genomes_Complete/${i}.fna.gz RefSeq/Genomes_Complete/temp_${i}.fna.gz \; \
    zcat RefSeq/Genomes_Complete/temp_${i}.fna.gz \| \
    sed \"s/^\>/\>${i}__/ \; s/ .*//\" \| \
    gzip \> RefSeq/Genomes_Complete/${i}.fna.gz \; \
    rm -f RefSeq/Genomes_Complete/temp_${i}.fna.gz
done > RefSeq/do_edit.sh
parallel -a RefSeq/do_edit.sh -j ${VAR_PARALLEL_EDIT}
rm -f RefSeq/do_edit.sh

for i in $(find RefSeq/Genomes_Incomplete/ -type f -name *.fna.gz | sed 's/RefSeq\/Genomes_Incomplete\///' | sed 's/.fna.gz//')
do
    echo mv RefSeq/Genomes_Incomplete/${i}.fna.gz RefSeq/Genomes_Incomplete/temp_${i}.fna.gz \; \
    zcat RefSeq/Genomes_Incomplete/temp_${i}.fna.gz \| \
    sed \"s/^\>/\>${i}__/ \; s/ .*//\" \| \
    gzip \> RefSeq/Genomes_Incomplete/${i}.fna.gz \; \
    rm -f RefSeq/Genomes_Incomplete/temp_${i}.fna.gz
done > RefSeq/do_edit.sh
parallel -a RefSeq/do_edit.sh -j ${VAR_PARALLEL_EDIT}
rm -f RefSeq/do_edit.sh

# make directories to store annotation information
mkdir RefSeq/Annotation_Complete/
mkdir RefSeq/Annotation_Incomplete/

# download annotation
awk 'BEGIN{FS=OFS="/";filesuffix="genomic.gff.gz"}{ftpdir=$0;asm=$10;file=asm"_"filesuffix;print "wget -P RefSeq/Annotation_Complete/ "ftpdir,file}' \
    RefSeq/ftpdirpathscomplete > RefSeq/ftpdownload
parallel -a RefSeq/ftpdownload -j ${VAR_PARALLEL_DL}
rm -f RefSeq/ftpdownload

awk 'BEGIN{FS=OFS="/";filesuffix="genomic.gff.gz"}{ftpdir=$0;asm=$10;file=asm"_"filesuffix;print "wget -P RefSeq/Annotation_Incomplete/ "ftpdir,file}' \
    RefSeq/ftpdirpathsincomplete > RefSeq/ftpdownload
parallel -a RefSeq/ftpdownload -j ${VAR_PARALLEL_DL}
rm -f RefSeq/ftpdownload

rm -f RefSeq/ftpdirpathscomplete
rm -f RefSeq/ftpdirpathsincomplete

# simplify annotation filenames
find RefSeq/Annotation_Complete/ -type f -name *.gff.gz | \
    awk -F '.[0-9]_' '{print "mv "$0" "$1".gff.gz"}' | bash
find RefSeq/Annotation_Incomplete/ -type f -name *.gff.gz | \
    awk -F '.[0-9]_' '{print "mv "$0" "$1".gff.gz"}' | bash

# extract rRNA genes and edit annotation seqids to contain RefSeq accession and match fasta headers
for i in $(find RefSeq/Annotation_Complete/ -type f -name '*.gff.gz' | sed 's/RefSeq\/Annotation_Complete\/// ; s/\.gff\.gz//')
do
    echo zcat RefSeq/Annotation_Complete/${i}.gff.gz \| \
    awk \'\$3==\"rRNA\"\' \| \
    sed \"s/^/${i}__/\" \> RefSeq/Annotation_Complete/${i}_rrna.gff
done > RefSeq/getRNA.sh
parallel -a RefSeq/getRNA.sh -j ${VAR_PARALLEL_EDIT}
rm -f RefSeq/getRNA.sh

for i in $(find RefSeq/Annotation_Incomplete/ -type f -name '*.gff.gz' | sed 's/RefSeq\/Annotation_Incomplete\/// ; s/\.gff\.gz//')
do
    echo zcat RefSeq/Annotation_Incomplete/${i}.gff.gz \| \
    awk \'\$3==\"rRNA\"\' \| \
    sed \"s/^/${i}__/\" \> RefSeq/Annotation_Incomplete/${i}_rrna.gff
done > RefSeq/getRNA.sh
parallel -a RefSeq/getRNA.sh -j ${VAR_PARALLEL_EDIT}
rm -f RefSeq/getRNA.sh

# make directories to store outputs
mkdir RefSeq/Outputs_Complete/
mkdir RefSeq/Outputs_Incomplete/

# make single multifasta files for all assemblies and annotations
find RefSeq/Genomes_Complete/ -type f -name "*.fna.gz" \
    -exec zcat {} + > RefSeq/Outputs_Complete/combined.fna
find RefSeq/Genomes_Incomplete/ -type f -name "*.fna.gz" \
    -exec zcat {} + > RefSeq/Outputs_Incomplete/combined.fna
find RefSeq/Annotation_Complete/ -type f -name "*_rrna.gff" \
    -exec cat {} + > RefSeq/Outputs_Complete/combined_rrna.gff
find RefSeq/Annotation_Incomplete/ -type f -name "*_rrna.gff" \
    -exec cat {} + > RefSeq/Outputs_Incomplete/combined_rrna.gff

# calculate the length of all sequences in mulitfasta
# this info will be used later to discard operons which cross the start/end boundary of the chromosome
awk '/^>/ {if (seqlen) print seqlen;print;seqlen=0;next} {seqlen+=length($0)}END{print seqlen}' RefSeq/Outputs_Complete/combined.fna | \
    sed 's/\.[0-9] .*//' | \
    sed 's/^>//' | \
    paste - - > RefSeq/Outputs_Complete/seq_length.tsv
awk '/^>/ {if (seqlen) print seqlen;print;seqlen=0;next} {seqlen+=length($0)}END{print seqlen}' RefSeq/Outputs_Incomplete/combined.fna | \
    sed 's/\.[0-9] .*//' | \
    sed 's/^>//' | \
    paste - - > RefSeq/Outputs_Incomplete/seq_length.tsv

# run custom R script to filter and reformat annotation
# takes barrnap output as input
# for every 16S gene, identifies the closest 23S gene with the same seqid and strand values
# identifies unlined rRNA operons and those that cross the chromosome start/end boundary, removes these and writes their details to 'operons_unlinked_or_crossingboundary.tsv'
# writes coordinates of ITS regions to 'combined_ITS.gff'
# writes coordinates of 16S-ITS-23S regions full_operons.gff' for extraction later with bedtools
# writes 'operon_identifiers.txt' with the unique operon identifiers that are substituted into the fasta headers later
Rscript Scripts/refseq_complete_filter_make_gff.R
Rscript Scripts/refseq_incomplete_filter_make_gff.R

# create a multifasta file containing all rRNA operons
bedtools getfasta \
    -fi RefSeq/Outputs_Complete/combined.fna \
    -bed RefSeq/Outputs_Complete/full_operons.gff \
    -fo RefSeq/Outputs_Complete/full_operons.fna \
    -s
rm -f RefSeq/Outputs_Complete/combined.fna.fai

bedtools getfasta \
    -fi RefSeq/Outputs_Incomplete/combined.fna \
    -bed RefSeq/Outputs_Incomplete/full_operons.gff \
    -fo RefSeq/Outputs_Incomplete/full_operons.fna \
    -s
rm -f RefSeq/Outputs_Incomplete/combined.fna.fai

# remove combined assembly and annotation files to save space
# temporarily keeping individual assembly files to do pairwise ANI
rm -f RefSeq/Outputs_Complete/combined.fna RefSeq/Outputs_Incomplete/combined.fna
rm -f RefSeq/Outputs_Complete/combined_rrna.gff RefSeq/Outputs_Combined/combined_rrna.gff

# renames sequences in operon fasta file by giving each operon a unique identifier
cut -f 2 RefSeq/Outputs_Complete/operon_identifiers.txt | \
    sed 's/^/>/' > RefSeq/Outputs_Complete/operon_identifiers_headers.txt
cut -f 2 RefSeq/Outputs_Incomplete/operon_identifiers.txt | \
    sed 's/^/>/' > RefSeq/Outputs_Incomplete/operon_identifiers_headers.txt

# replaces the fasta headers with unique operon identifiers
sed -i 's/>.* />/' RefSeq/Outputs_Complete/operon_identifiers_headers.txt
awk 'NR%2==0' RefSeq/Outputs_Complete/full_operons.fna | \
    paste -d '\n' RefSeq/Outputs_Complete/operon_identifiers_headers.txt - > RefSeq/Outputs_Complete/full_operons_identifiers.fna

sed -i 's/>.* />/' RefSeq/Outputs_Incomplete/operon_identifiers_headers.txt
awk 'NR%2==0' RefSeq/Outputs_Incomplete/full_operons.fna | \
    paste -d '\n' RefSeq/Outputs_Incomplete/operon_identifiers_headers.txt - > RefSeq/Outputs_Incomplete/full_operons_identifiers.fna

# create list of assembly_accessions and taxonIDs from summary file
zcat RefSeq/prokaryote_assembly_summary_latest.txt.gz | \
    cut -f 1,6 > RefSeq/assemblyaccession_taxid.txt

# modify first column to remove the .[0-9] at the end of accession to match prefix added to fasta header and seqid
sed -i 's/\.[0-9]//' RefSeq/assemblyaccession_taxid.txt

# assigning taxonomy lineage to taxids
# output is a tab delim file with one operon per row
taxonkit lineage \
    --data-dir ${VAR_DB_TAXONDB} RefSeq/assemblyaccession_taxid.txt \
    -i 2 | \
    awk '$2!=""' > RefSeq/assemblyaccession_taxid_lineage.txt

# converting operon lineage info to MPA-style
taxonkit reformat \
    -a \
    -i 3 \
    -d ";" \
    -f "{k}|{p}|{c}|{o}|{f}|{g}|{s}|{t}" \
    -P RefSeq/assemblyaccession_taxid_lineage.txt \
    --data-dir ${VAR_DB_TAXONDB} | \
    cut -f 1,4 > RefSeq/assemblyaccession_taxid_lineage_mpa.txt
sed -i 's/^>//' RefSeq/Outputs_Complete/operon_identifier_lineage_mpa.txt

# tidying up intermediary files
rm -f RefSeq/assemblyaccession_taxid.txt RefSeq/assemblyaccession_taxid_lineage.txt
rm -f RefSeq/Outputs_Complete/operon_identifiers_headers.txt  RefSeq/Outputs_Complete/full_operons.fna
rm -f RefSeq/Outputs_Incomplete/operon_identifiers_headers.txt RefSeq/Outputs_Incomplete/full_operons.fna

# sorts sequences by length (important for NR clustering, will retain longest sequence in each cluster)
sortbyname.sh \
    in=RefSeq/Outputs_Complete/full_operons_identifiers.fna \
    out=RefSeq/Outputs_Complete/full_operons_sorted.fna \
    length \
    descending
rm -f RefSeq/Outputs_Complete/full_operons_identifiers.fna

sortbyname.sh \
    in=RefSeq/Outputs_Incomplete/full_operons_identifiers.fna \
    out=RefSeq/Outputs_Incomplete/full_operons_sorted.fna \
    length \
    descending
rm -f RefSeq/Outputs_Incomplete/full_operons_identifiers.fna

# removes redundant sequences from Complete Assemblies dataset based on specified identity cutoff
# also extracts relevant cluster information
vsearch \
    --cluster_smallmem RefSeq/Outputs_Complete/full_operons_sorted.fna \
    --id ${VAR_CLUSTERID_VSEARCH} \
    --centroids RefSeq/Outputs_Complete/nrRep.fna \
    --uc RefSeq/Outputs_Complete/clusters.tsv \
    --threads ${VAR_THREADS_VSEARCH} \
    -log RefSeq/Outputs_Complete/vsearchlog.txt \
    --consout RefSeq/Outputs_Complete/nrCon.fna
cat RefSeq/Outputs_Complete/clusters.tsv | \
    awk '$1=="S"' > RefSeq/Outputs_Complete/vsearch_centroids.tsv
cat RefSeq/Outputs_Complete/clusters.tsv | \
    awk '$1=="C"' > RefSeq/Outputs_Complete/vsearch_clusters.tsv
cat RefSeq/Outputs_Complete/clusters.tsv | \
    awk '$1=="H"' > RefSeq/Outputs_Complete/vsearch_hits.tsv

# adds operon sequences from Incomplete Assemblies to the NR dataset generated from complete assemblies, reruns clustering, and extracts relevant cluster information
# only adds NR sequences from Incomplete genome assemblies if they are not already represented in the Complete dataset
# also combines the taxonomy tables of these incomplete operons
mkdir RefSeq/Outputs_Combined/
cat RefSeq/Outputs_Complete/nrRep.fna RefSeq/Outputs_Incomplete/full_operons_sorted.fna > RefSeq/Outputs_Combined/full_operons_sorted.fna
vsearch \
    --cluster_smallmem RefSeq/Outputs_Combined/full_operons_sorted.fna \
    --id ${VAR_CLUSTERID_VSEARCH} \
    --centroids RefSeq/Outputs_Combined/nrRep.fna \
    --uc RefSeq/Outputs_Combined/clusters.tsv \
    --threads ${VAR_THREADS_VSEARCH} \
    -log RefSeq/Outputs_Combined/vsearchlog.txt \
    --consout RefSeq/Outputs_Combined/nrCon.fna \
    --usersort
cat RefSeq/Outputs_Combined/clusters.tsv | \
    awk '$1=="S"' > RefSeq/Outputs_Combined/vsearch_centroids.tsv
cat RefSeq/Outputs_Combined/clusters.tsv | \
    awk '$1=="C"' > RefSeq/Outputs_Combined/vsearch_clusters.tsv
cat RefSeq/Outputs_Combined/clusters.tsv | \
    awk '$1=="H"' > RefSeq/Outputs_Combined/vsearch_hits.tsv

# edit fasta headers of nrCon sequence files
sed -i 's/centroid=// ; s/;seqs=.*//' RefSeq/Outputs_Complete/nrCon.fna 
sed -i 's/centroid=// ; s/;seqs=.*//' RefSeq/Outputs_Combined/nrCon.fna

# move database files to single directory
mkdir -p Database/Complete/
mkdir -p Database/Combined/

sortbyname.sh \
    in=RefSeq/Outputs_Complete/full_operons_sorted.fna \
    out=Database/Complete/full.fna.gz
sortbyname.sh \
    in=RefSeq/Outputs_Combined/full_operons_sorted.fna \
    out=Database/Combined/full.fna.gz

sortbyname.sh \
    in=RefSeq/Outputs_Complete/nrRep.fna \
    out=Database/Complete/nrRep.fna.gz
sortbyname.sh \
    in=RefSeq/Outputs_Complete/nrCon.fna \
    out=Database/Complete/nrCon.fna.gz
sortbyname.sh \
    in=RefSeq/Outputs_Combined/nrRep.fna \
    out=Database/Combined/nrRep.fna.gz
sortbyname.sh \
    in=RefSeq/Outputs_Combined/nrCon.fna \
    out=Database/Combined/nrCon.fna.gz

# remove database sequence files from output directories
rm -f RefSeq/Outputs_Complete/full_operons_sorted.fna
rm -f RefSeq/Outputs_Combined/full_operons_sorted.fna
rm -f RefSeq/Outputs_Complete/nrRep.fna
rm -f RefSeq/Outputs_Complete/nrCon.fna
rm -f RefSeq/Outputs_Combined/nrRep.fna
rm -f RefSeq/Outputs_Combined/nrCon.fna

# generate database taxonomy files and move to database directory
Rscript Scripts/refseq_get_taxonomy.R
mv Outputs_Complete/taxFull.tsv Database/Complete/taxFull.tsv
mv Outputs_Combined/taxFull.tsv Database/Combined/taxFull.tsv
mv Outputs_Complete/taxRep.tsv Database/Complete/taxRep.tsv
mv Outputs_Complete/taxLCA.tsv Database/Complete/taxLCA.tsv
mv Outputs_Complete/taxMaj.tsv Database/Complete/taxMaj.tsv
mv Outputs_Combined/taxRep.tsv Database/Combined/taxRep.tsv
mv Outputs_Combined/taxLCA.tsv Database/Combined/taxLCA.tsv
mv Outputs_Combined/taxMaj.tsv Database/Combined/taxMaj.tsv

