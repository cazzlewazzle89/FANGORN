#!/bin/sh
VAR_PARALLEL_DL=50 # number of parallel processes to run when downloading assemblies & annotations
VAR_PARALLEL_EDIT=50 # number of parallel processes to run when editing seqid headers
VAR_THREADS_VSEARCH=50 # number of threads to use then clustering operon sequences with vsearch
VAR_CLUSTERID_VSEARCH=0.999 # vsearch clustering identity
VAR_BARRNAP_THREADS=20
VAR_OUTPUT_DIRECTORY=${PWD}/Fangorn_GTDB # specify output directory name for databases

# make output directory and move shell to run from there
mkdir ${VAR_OUTPUT_DIRECTORY}
cd ${VAR_OUTPUT_DIRECTORY}

# make directories to store genomes used for building the database
mkdir -p Genomes_Complete/
mkdir -p Genomes_Incomplete/

# download archaea and bacteria taxonomy and combine into a single file
wget https://data.gtdb.ecogenomic.org/releases/latest/ar122_taxonomy.tsv.gz \
    -O GTDB/ar122_taxonomy.tsv.gz \
    --no-check-certificate
gunzip GTDB/ar122_taxonomy.tsv.gz

wget https://data.gtdb.ecogenomic.org/releases/latest/bac120_taxonomy.tsv.gz \
    -O GTDB/bac120_taxonomy.tsv.gz \
    --no-check-certificate
gunzip GTDB/bac120_taxonomy.tsv.gz

mv GTDB/ar122_taxonomy.tsv GTDB/gtdb_taxonomy.tsv
sed -e '1d' GTDB/bac120_taxonomy.tsv >> GTDB/gtdb_taxonomy.tsv
rm -f GTDB/bac120_taxonomy.tsv

# download archaea and bacteria metadata and combine into a single file
wget https://data.gtdb.ecogenomic.org/releases/latest/ar122_metadata.tar.gz \
    -O GTDB/ar122_metadata.tar.gz \
    --no-check-certificate
tar -xzvf GTDB/ar122_metadata.tar.gz -C GTDB/
rm -f GTDB/ar122_metadata.tar.gz

wget https://data.gtdb.ecogenomic.org/releases/latest/bac120_metadata.tar.gz \
    -O GTDB/bac120_metadata.tar.gz \
    --no-check-certificate
tar -xzvf GTDB/bac120_metadata.tar.gz  -C GTDB/
rm -f GTDB/bac120_metadata.tar.gz

mv GTDB/ar122_metadata_*.tsv GTDB/gtdb_metadata.tsv
sed -e '1d' GTDB/bac120_metadata_*.tsv >> GTDB/gtdb_metadata.tsv
rm -f GTDB/bac120_metadata_*.tsv

# download genbank and refseq assembly information
wget ftp://ftp.ncbi.nlm.nih.gov/genomes/ASSEMBLY_REPORTS/assembly_summary_genbank.txt \
    -O GTDB/assembly_summary_genbank.txt
wget ftp://ftp.ncbi.nlm.nih.gov/genomes/ASSEMBLY_REPORTS/assembly_summary_refseq.txt \
    -O GTDB/assembly_summary_refseq.txt

# get assembly info for each assembly in the GTDB database
grep '^GB_GCA' GTDB/gtdb_taxonomy.tsv | \
    cut -f 1 | \
    sed 's/GB_//' > GTDB/genbank_gtdb_accessions.txt
grep '^RS_GCF' GTDB/gtdb_taxonomy.tsv \
    | cut -f 1 | \
    sed 's/RS_//' > GTDB/refseq_gtdb_accessions.txt
Rscript Scripts/gtdb_get_assembly_info.R

# extract ftp directory paths for each assembly
cut -f 20 GTDB/assembly_summary_genbank_gtdb_complete.txt > GTDB/ftpdirpaths_genbank_complete
cut -f 20 GTDB/assembly_summary_genbank_gtdb_incomplete.txt > GTDB/ftpdirpaths_genbank_incomplete
cut -f 20 GTDB/assembly_summary_refseq_gtdb_complete.txt > GTDB/ftpdirpaths_refseq_complete
cut -f 20 GTDB/assembly_summary_refseq_gtdb_incomplete.txt > GTDB/ftpdirpaths_refseq_incomplete

# compress assembly information
# retaining to keep information about database version
gzip GTDB/assembly_summary_genbank.txt 
gzip GTDB/assembly_summary_refseq.txt 
gzip GTDB/assembly_summary_genbank_gtdb_complete.txt 
gzip GTDB/assembly_summary_genbank_gtdb_incomplete.txt 
gzip GTDB/assembly_summary_refseq_gtdb_complete.txt
gzip GTDB/assembly_summary_refseq_gtdb_incomplete.txt

# remove files listing assession numbers
rm -f GTDB/genbank_gtdb_accessions.txt GTDB/refseq_gtdb_accessions.txt

# make directories to store downloaded genomes
mkdir GTDB/Genomes_Complete/
mkdir GTDB/Genomes_Incomplete/

# download nucleotide sequences
awk 'BEGIN{FS=OFS="/";filesuffix="genomic.fna.gz"}{ftpdir=$0;asm=$10;file=asm"_"filesuffix;print "wget -P GTDB/Genomes_Complete/ "ftpdir,file}' \
    GTDB/ftpdirpaths_genbank_complete > GTDB/ftpdownload
parallel -a GTDB/ftpdownload -j ${VAR_PARALLEL_DL}
rm -f GTDB/ftpdownload

awk 'BEGIN{FS=OFS="/";filesuffix="genomic.fna.gz"}{ftpdir=$0;asm=$10;file=asm"_"filesuffix;print "wget -P GTDB/Genomes_Complete/ "ftpdir,file}' \
    GTDB/ftpdirpaths_refseq_complete > GTDB/ftpdownload
parallel -a GTDB/ftpdownload -j ${VAR_PARALLEL_DL}
rm -f GTDB/ftpdownload

awk 'BEGIN{FS=OFS="/";filesuffix="genomic.fna.gz"}{ftpdir=$0;asm=$10;file=asm"_"filesuffix;print "wget -P GTDB/Genomes_Incomplete/ "ftpdir,file}' \
    GTDB/ftpdirpaths_genbank_incomplete > GTDB/ftpdownload
parallel -a GTDB/ftpdownload -j ${VAR_PARALLEL_DL}
rm -f GTDB/ftpdownload

awk 'BEGIN{FS=OFS="/";filesuffix="genomic.fna.gz"}{ftpdir=$0;asm=$10;file=asm"_"filesuffix;print "wget -P GTDB/Genomes_Incomplete/ "ftpdir,file}' \
    GTDB/ftpdirpaths_refseq_incomplete > GTDB/ftpdownload
parallel -a GTDB/ftpdownload -j ${VAR_PARALLEL_DL}
rm -f GTDB/ftpdownload

# simplify assembly filenames
find GTDB/Genomes_Complete/ -type f -name '*.fna.gz' | \
    awk -F '.[0-9]_' '{print "mv "$0" "$1".fna.gz"}' > rename.sh
sh rename.sh

find GTDB/Genomes_Incomplete/ -type f -name '*.fna.gz' | \
    awk -F '.[0-9]_' '{print "mv "$0" "$1".fna.gz"}' > rename.sh
sh rename.sh
rm -f rename.sh

# edit fasta headers to contain the RefSeq or GenBank assembly accession
find GTDB/Genomes_Complete/ -type f -name '*.fna.gz' | sed 's/GTDB\/Genomes_Complete\///' | sed 's/.fna.gz//' > GTDB/temp_filelist
for i in $(cat GTDB/temp_filelist)
do
    echo mv GTDB/Genomes_Complete/${i}.fna.gz GTDB/Genomes_Complete/temp_${i}.fna.gz \; \
    zcat GTDB/Genomes_Complete/temp_${i}.fna.gz \| \
    sed \"s/^\>/\>${i}__/ \; s/ .*//\" \| \
    gzip \> GTDB/Genomes_Complete/${i}.fna.gz \; \
    rm -f GTDB/Genomes_Complete/temp_${i}.fna.gz
done > GTDB/do_edit.sh
parallel -a GTDB/do_edit.sh -j ${VAR_PARALLEL_EDIT}
rm -f GTDB/do_edit.sh GTDB/temp_filelist

find GTDB/Genomes_Incomplete/ -type f -name '*.fna.gz' | sed 's/GTDB\/Genomes_Incomplete\///' | sed 's/\.fna\.gz//' > GTDB/temp_filelist
for i in $(cat GTDB/temp_filelist)
do
    echo mv GTDB/Genomes_Incomplete/${i}.fna.gz GTDB/Genomes_Incomplete/temp_${i}.fna.gz \; \
    zcat GTDB/Genomes_Incomplete/temp_${i}.fna.gz \| \
    sed \"s/^\>/\>${i}__/ \; s/ .*//\" \| \
    gzip \> GTDB/Genomes_Incomplete/${i}.fna.gz \; \
    rm -f GTDB/Genomes_Incomplete/temp_${i}.fna.gz
done > GTDB/do_edit.sh
parallel -a GTDB/do_edit.sh -j ${VAR_PARALLEL}
rm -f GTDB/do_edit.sh GTDB/temp_filelist

# make directories to store downloaded annotation information
mkdir GTDB/Annotation_Complete
mkdir GTDB/Annotation_Incomplete

# download annotation information
awk 'BEGIN{FS=OFS="/";filesuffix="genomic.gff.gz"}{ftpdir=$0;asm=$10;file=asm"_"filesuffix;print "wget -P GTDB/Annotation_Complete/ "ftpdir,file}' \
    GTDB/ftpdirpaths_genbank_complete > GTDB/ftpdownload
parallel -a GTDB/ftpdownload -j ${VAR_PARALLEL_DL}
rm -f GTDB/ftpdownload

awk 'BEGIN{FS=OFS="/";filesuffix="genomic.gff.gz"}{ftpdir=$0;asm=$10;file=asm"_"filesuffix;print "wget -P GTDB/Annotation_Complete/ "ftpdir,file}' \
    GTDB/ftpdirpaths_refseq_complete > GTDB/ftpdownload
parallel -a GTDB/ftpdownload -j ${VAR_PARALLEL_DL}
rm -f GTDB/ftpdownload

awk 'BEGIN{FS=OFS="/";filesuffix="genomic.gff.gz"}{ftpdir=$0;asm=$10;file=asm"_"filesuffix;print "wget -P GTDB/Annotation_Incomplete/ "ftpdir,file}' \
    GTDB/ftpdirpaths_genbank_incomplete > GTDB/ftpdownload
parallel -a GTDB/ftpdownload -j ${VAR_PARALLEL_DL}
rm -f GTDB/ftpdownload

awk 'BEGIN{FS=OFS="/";filesuffix="genomic.gff.gz"}{ftpdir=$0;asm=$10;file=asm"_"filesuffix;print "wget -P GTDB/Annotation_Incomplete/ "ftpdir,file}' \
GTDB/ftpdirpaths_refseq_incomplete > GTDB/ftpdownload
parallel -a GTDB/ftpdownload -j ${VAR_PARALLEL_DL}
rm -f GTDB/ftpdownload

# simplify annotation filenames
find GTDB/Annotation_Complete/ -type f -name '*.gff.gz' | \
    awk -F '.[0-9]_' '{print "mv "$0" "$1".gff.gz"}' > rename.sh
sh rename.sh
find GTDB/Annotation_Incomplete/ -type f -name '*.gff.gz' | \
    awk -F '.[0-9]_' '{print "mv "$0" "$1".gff.gz"}' > rename.sh
sh rename.sh
rm -f rename.sh

# identify unannotated genomes
Rscript Scripts/genomes_missing_annotation.R

# identify rRNA genes in unannotated genomes
for i in $(cat GTDB/missingannotation_complete.txt)
do
    zcat GTDB/Genomes_Complete/${i}.fna.gz | \
    barrnap --threads ${VAR_BARRNAP_THREADS} \
        --reject 0.8 > GTDB/Annotation_Complete/${i}_rrna.gff
done

for i in $(cat GTDB/missingannotation_incomplete.txt)
do
    zcat GTDB/Genomes_Incomplete/${i}.fna.gz | \
        barrnap --threads ${VAR_BARRNAP_THREADS} \
            --reject 0.8 > GTDB/Annotation_Incomplete/${i}_rrna.gff
done

# extract rRNA genes and edit annotation seqids to contain RefSeq or GenBank accession and match fasta headers
for i in $(find GTDB/Annotation_Complete/ -type f -name '*.gff.gz' | sed 's/GTDB\/Annotation_Complete\/// ; s/\.gff\.gz//')
do
    echo zcat GTDB/Annotation_Complete/${i}.gff.gz \| \
    awk \'\$3==\"rRNA\"\' \| \
    sed \"s/^/${i}__/\" \> GTDB/Annotation_Complete/${i}_rrna.gff
done > GTDB/getRNA.sh
parallel -a GTDB/getRNA.sh -j ${VAR_PARALLEL_EDIT}
rm -f GTDB/getRNA.sh

for i in $(find GTDB/Annotation_Incomplete/ -type f -name '*.gff.gz' | sed 's/GTDB\/Annotation_Incomplete\/// ; s/\.gff\.gz//')
do
    echo zcat GTDB/Annotation_Incomplete/${i}.gff.gz \| \
    awk \'\$3==\"rRNA\"\' \| \
    sed \"s/^/${i}__/\" \> \
    GTDB/Annotation_Incomplete/${i}_rrna.gff
done > GTDB/getRNA.sh
parallel -a GTDB/getRNA.sh -j ${VAR_PARALLEL_EDIT}
rm -f GTDB/getRNA.sh

# make directories to store outputs
mkdir GTDB/Outputs_Complete/
mkdir GTDB/Outputs_Incomplete/

# make single multifasta files for all assemblies and annotations
find GTDB/Genomes_Complete/ -type f -name "*.fna.gz" \
    -exec zcat {} + > GTDB/Outputs_Complete/combined.fna
find GTDB/Genomes_Incomplete/ -type f -name "*.fna.gz" \
    -exec zcat {} + > GTDB/Outputs_Incomplete/combined.fna
find GTDB/Annotation_Complete/ -type f -name "*_rrna.gff" \
    -exec cat {} + > GTDB/Outputs_Complete/combined_rrna.gff
find GTDB/Annotation_Incomplete/ -type f -name "*_rrna.gff" \
    -exec cat {} + > GTDB/Outputs_Incomplete/combined_rrna.gff

# calculate the length of all sequences in mulitfasta
# this info will be used later to discard operons which cross the start/end boundary of the chromosome
awk '/^>/ {if (seqlen) print seqlen;print;seqlen=0;next} {seqlen+=length($0)}END{print seqlen}' \
    GTDB/Outputs_Complete/combined.fna | \
    sed 's/\.[0-9] .*//' | \
    sed 's/^>//' | \
    paste - - > GTDB/Outputs_Complete/seq_length.tsv

awk '/^>/ {if (seqlen) print seqlen;print;seqlen=0;next} {seqlen+=length($0)}END{print seqlen}' \
    GTDB/Outputs_Incomplete/combined.fna | \
    sed 's/\.[0-9] .*//' | \
    sed 's/^>//' | \
    paste - - > GTDB/Outputs_Incomplete/seq_length.tsv

# run custom R script to filter and reformat annotation
# takes barrnap output as input
# for every 16S gene, identifies the closest 23S gene with the same seqid and strand values
# identifies unlined rRNA operons and those that cross the chromosome start/end boundary, removes these and writes their details to 'operons_unlinked_or_crossingboundary.tsv'
# writes coordinates of ITS regions to 'combined_ITS.gff'
# writes coordinates of 16S-ITS-23S regions full_operons.gff' for extraction later with bedtools
# writes 'operon_identifiers.txt' with the unique operon identifiers that are substituted into the fasta headers later
Rscript Scripts/gtdb_complete_filter_make_gff.R
Rscript Scripts/gtdb_incomplete_filter_make_gff.R.R

# create a multifasta file containing all rRNA operons
bedtools getfasta \
    -fi GTDB/Outputs_Complete/combined.fna \
    -bed GTDB/Outputs_Complete/full_operons.gff \
    -fo GTDB/Outputs_Complete/full_operons.fna \
    -s
rm -f GTDB/Outputs_Complete/combined.fna.fai

bedtools getfasta \
    -fi GTDB/Outputs_Incomplete/combined.fna \
    -bed GTDB/Outputs_Incomplete/full_operons.gff \
    -fo GTDB/Outputs_Incomplete/full_operons.fna \
    -s
rm -f GTDB/Outputs_Incomplete/combined.fna.fai

# remove combined assembly and annotation files to save space
# temporarily keeping individual assembly files to do pairwise ANI
rm -f GTDB/Outputs_Complete/combined.fna GTDB/Outputs_Incomplete/combined.fna
rm -f GTDB/Outputs_Complete/combined_rrna.gff GTDB/Outputs_Combined/combined_rrna.gff

# renames sequences in operon fasta file by giving each operon a unique identifier
cut -f 2 GTDB/Outputs_Complete/operon_identifiers.txt | \
    sed 's/^/>/' > GTDB/Outputs_Complete/operon_identifiers_headers.txt
cut -f 2 GTDB/Outputs_Incomplete/operon_identifiers.txt | \
    sed 's/^/>/' > GTDB/Outputs_Incomplete/operon_identifiers_headers.txt

# replaces the fasta headers with unique operon identifiers
sed -i 's/>.* />/' GTDB/Outputs_Complete/operon_identifiers_headers.txt
awk 'NR%2==0' GTDB/Outputs_Complete/full_operons.fna | \
    paste -d '\n' GTDB/Outputs_Complete/operon_identifiers_headers.txt - > GTDB/Outputs_Complete/full_operons_identifiers.fna
sed -i 's/>.* />/' GTDB/Outputs_Incomplete/operon_identifiers_headers.txt
awk 'NR%2==0' GTDB/Outputs_Incomplete/full_operons.fna | \
    paste -d '\n' GTDB/Outputs_Incomplete/operon_identifiers_headers.txt - > RefSeq/Outputs_Incomplete/full_operons_identifiers.fna

# sorts sequences by length (important for NR clustering, will retain longest sequence in each cluster)
sortbyname.sh \
    in=GTDB/Outputs_Complete/full_operons_identifiers.fna \
    out=GTDB/Outputs_Complete/full_operons_sorted.fna \
    length \
    descending
rm -f GTDB/Outputs_Complete/full_operons_identifiers.fna

sortbyname.sh \
    in=GTDB/Outputs_Incomplete/full_operons_identifiers.fna \
    out=GTDB/Outputs_Incomplete/full_operons_sorted.fna \
    length \
    descending
rm -f GTDB/Outputs_Incomplete/full_operons_identifiers.fna

# removes redundant sequences from Complete Assemblies dataset based on specified identity cutoff
# also extracts relevant cluster information
vsearch \
    --cluster_smallmem GTDB/Outputs_Complete/full_operons_sorted.fna \
    --id ${VAR_CLUSTERID_VSEARCH} \
    --centroids GTDB/Outputs_Complete/nrRep.fna \
    --uc GTDB/Outputs_Complete/clusters.tsv \
    --threads ${VAR_THREADS_VSEARCH} \
    -log GTDB/Outputs_Complete/vsearchlog.txt \
    --consout GTDB/Outputs_Complete/nrCon.fna
cat GTDB/Outputs_Complete/clusters.tsv | \
    awk '$1=="S"' > GTDB/Outputs_Complete/vsearch_centroids.tsv
cat GTDB/Outputs_Complete/clusters.tsv | \
    awk '$1=="C"' > GTDB/Outputs_Complete/vsearch_clusters.tsv
cat GTDB/Outputs_Complete/clusters.tsv | \
    awk '$1=="H"' > GTDB/Outputs_Complete/vsearch_hits.tsv

# adds operon sequences from Incomplete Assemblies to the NR dataset generated from complete assemblies, reruns clustering, and extracts relevant cluster information
# only adds NR sequences from Incomplete genome assemblies if they are not already represented in the Complete dataset
# also combines the taxonomy tables of these incomplete operons
mkdir GTDB/Outputs_Combined/
cat GTDB/Outputs_Complete/nrRep.fna GTDB/Outputs_Incomplete/full_operons_sorted.fna > GTDB/Outputs_Combined/full_operons_sorted.fna
vsearch \
    --cluster_smallmem GTDB/Outputs_Combined/full_operons_sorted.fna \
    --id ${VAR_CLUSTERID_VSEARCH} \
    --centroids GTDB/Outputs_Combined/nrRep.fna \
    --uc GTDB/Outputs_Combined/clusters.tsv \
    --threads ${VAR_THREADS_VSEARCH} \
    -log GTDB/Outputs_Combined/vsearchlog.txt \
    --consout GTDB/Outputs_Combined/nrCon.fna \
    --usersort
cat GTDB/Outputs_Combined/clusters.tsv | \
    awk '$1=="S"' > GTDB/Outputs_Combined/vsearch_centroids.tsv
cat GTDB/Outputs_Combined/clusters.tsv | \
    awk '$1=="C"' > GTDB/Outputs_Combined/vsearch_clusters.tsv
cat GTDB/Outputs_Combined/clusters.tsv | \
    awk '$1=="H"' > GTDB/Outputs_Combined/vsearch_hits.tsv

# edit fasta headers of nrCon sequence files
sed -i 's/centroid=// ; s/;seqs=.*//' GTDB/Outputs_Complete/nrCon.fna 
sed -i 's/centroid=// ; s/;seqs=.*//' GTDB/Outputs_Combined/nrCon.fna

# move database files to single directory
mkdir -p Database/Complete/
mkdir -p Database/Combined/

sortbyname.sh \
    in=GTDB/Outputs_Complete/full_operons_sorted.fna \
    out=Database/Complete/full.fna.gz
sortbyname.sh \
    in=GTDB/Outputs_Combined/full_operons_sorted.fna \
    out=Database/Combined/full.fna.gz

sortbyname.sh \
    in=GTDB/Outputs_Complete/nrRep.fna \
    out=Database/Complete/nrRep.fna.gz
sortbyname.sh \
    in=GTDB/Outputs_Complete/nrCon.fna \
    out=Database/Complete/nrCon.fna.gz
sortbyname.sh \
    in=GTDB/Outputs_Combined/nrRep.fna \
    out=Database/Combined/nrRep.fna.gz
sortbyname.sh \
    in=GTDB/Outputs_Combined/nrCon.fna \
    out=Database/Combined/nrCon.fna.gz

# remove database sequence files from output directories
rm -f GTDB/Outputs_Complete/full_operons_sorted.fna
rm -f GTDB/Outputs_Combined/full_operons_sorted.fna
rm -f GTDB/Outputs_Complete/nrRep.fna
rm -f GTDB/Outputs_Complete/nrCon.fna
rm -f GTDB/Outputs_Combined/nrRep.fna
rm -f GTDB/Outputs_Combined/nrCon.fna

# get GTDB taxonomy for each rrn operon in databases
Rscript Scripts/gtdb_get_taxonomy.R
mv GTDB/Outputs_Complete/taxFull.tsv Database/Complete/taxFull.tsv
mv GTDB/Outputs_Combined/taxFull.tsv Database/Combined/taxFull.tsv
mv GTDB/Outputs_Complete/taxRep.tsv Database/Complete/taxRep.tsv
mv GTDB/Outputs_Complete/taxLCA.tsv Database/Complete/taxLCA.tsv
mv GTDB/Outputs_Complete/taxMaj.tsv Database/Complete/taxMaj.tsv
mv GTDB/Outputs_Combined/taxRep.tsv Database/Combined/taxRep.tsv
mv GTDB/Outputs_Combined/taxLCA.tsv Database/Combined/taxLCA.tsv
mv GTDB/Outputs_Combined/taxMaj.tsv Database/Combined/taxMaj.tsv

