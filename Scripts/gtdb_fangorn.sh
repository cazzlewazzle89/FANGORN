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
    -O ar122_taxonomy.tsv.gz \
    --no-check-certificate
gunzip ar122_taxonomy.tsv.gz

wget https://data.gtdb.ecogenomic.org/releases/latest/bac120_taxonomy.tsv.gz \
    -O bac120_taxonomy.tsv.gz \
    --no-check-certificate
gunzip bac120_taxonomy.tsv.gz

mv ar122_taxonomy.tsv gtdb_taxonomy.tsv
sed -e '1d' bac120_taxonomy.tsv >> gtdb_taxonomy.tsv
rm -f bac120_taxonomy.tsv

# download archaea and bacteria metadata and combine into a single file
wget https://data.gtdb.ecogenomic.org/releases/latest/ar122_metadata.tar.gz \
    -O ar122_metadata.tar.gz \
    --no-check-certificate
tar -xzvf ar122_metadata.tar.gz -C 
rm -f ar122_metadata.tar.gz

wget https://data.gtdb.ecogenomic.org/releases/latest/bac120_metadata.tar.gz \
    -O bac120_metadata.tar.gz \
    --no-check-certificate
tar -xzvf bac120_metadata.tar.gz  -C 
rm -f bac120_metadata.tar.gz

mv ar122_metadata_*.tsv gtdb_metadata.tsv
sed -e '1d' bac120_metadata_*.tsv >> gtdb_metadata.tsv
rm -f bac120_metadata_*.tsv

# download genbank and refseq assembly information
wget ftp://ftp.ncbi.nlm.nih.gov/genomes/ASSEMBLY_REPORTS/assembly_summary_genbank.txt \
    -O assembly_summary_genbank.txt
wget ftp://ftp.ncbi.nlm.nih.gov/genomes/ASSEMBLY_REPORTS/assembly_summary_refseq.txt \
    -O assembly_summary_refseq.txt

# get assembly info for each assembly in the GTDB database
grep '^GB_GCA' gtdb_taxonomy.tsv | \
    cut -f 1 | \
    sed 's/GB_//' > genbank_gtdb_accessions.txt
grep '^RS_GCF' gtdb_taxonomy.tsv \
    | cut -f 1 | \
    sed 's/RS_//' > refseq_gtdb_accessions.txt
Rscript Scripts/gtdb_get_assembly_info.R

# extract ftp directory paths for each assembly
cut -f 20 assembly_summary_genbank_gtdb_complete.txt > ftpdirpaths_genbank_complete
cut -f 20 assembly_summary_genbank_gtdb_incomplete.txt > ftpdirpaths_genbank_incomplete
cut -f 20 assembly_summary_refseq_gtdb_complete.txt > ftpdirpaths_refseq_complete
cut -f 20 assembly_summary_refseq_gtdb_incomplete.txt > ftpdirpaths_refseq_incomplete

# compress assembly information
# retaining to keep information about database version
gzip assembly_summary_genbank.txt 
gzip assembly_summary_refseq.txt 
gzip assembly_summary_genbank_gtdb_complete.txt 
gzip assembly_summary_genbank_gtdb_incomplete.txt 
gzip assembly_summary_refseq_gtdb_complete.txt
gzip assembly_summary_refseq_gtdb_incomplete.txt

# remove files listing assession numbers
rm -f genbank_gtdb_accessions.txt refseq_gtdb_accessions.txt

# make directories to store downloaded genomes
mkdir Genomes_Complete/
mkdir Genomes_Incomplete/

# download nucleotide sequences
awk 'BEGIN{FS=OFS="/";filesuffix="genomic.fna.gz"}{ftpdir=$0;asm=$10;file=asm"_"filesuffix;print "wget -P Genomes_Complete/ "ftpdir,file}' \
    ftpdirpaths_genbank_complete > ftpdownload
parallel -a ftpdownload -j ${VAR_PARALLEL_DL}
rm -f ftpdownload

awk 'BEGIN{FS=OFS="/";filesuffix="genomic.fna.gz"}{ftpdir=$0;asm=$10;file=asm"_"filesuffix;print "wget -P Genomes_Complete/ "ftpdir,file}' \
    ftpdirpaths_refseq_complete > ftpdownload
parallel -a ftpdownload -j ${VAR_PARALLEL_DL}
rm -f ftpdownload

awk 'BEGIN{FS=OFS="/";filesuffix="genomic.fna.gz"}{ftpdir=$0;asm=$10;file=asm"_"filesuffix;print "wget -P Genomes_Incomplete/ "ftpdir,file}' \
    ftpdirpaths_genbank_incomplete > ftpdownload
parallel -a ftpdownload -j ${VAR_PARALLEL_DL}
rm -f ftpdownload

awk 'BEGIN{FS=OFS="/";filesuffix="genomic.fna.gz"}{ftpdir=$0;asm=$10;file=asm"_"filesuffix;print "wget -P Genomes_Incomplete/ "ftpdir,file}' \
    ftpdirpaths_refseq_incomplete > ftpdownload
parallel -a ftpdownload -j ${VAR_PARALLEL_DL}
rm -f ftpdownload

# simplify assembly filenames
find Genomes_Complete/ -type f -name '*.fna.gz' | \
    awk -F '.[0-9]_' '{print "mv "$0" "$1".fna.gz"}' > rename.sh
sh rename.sh

find Genomes_Incomplete/ -type f -name '*.fna.gz' | \
    awk -F '.[0-9]_' '{print "mv "$0" "$1".fna.gz"}' > rename.sh
sh rename.sh
rm -f rename.sh

# edit fasta headers to contain the RefSeq or GenBank assembly accession
find Genomes_Complete/ -type f -name '*.fna.gz' | sed 's/GTDB\/Genomes_Complete\///' | sed 's/.fna.gz//' > temp_filelist
for i in $(cat temp_filelist)
do
    echo mv Genomes_Complete/${i}.fna.gz Genomes_Complete/temp_${i}.fna.gz \; \
    zcat Genomes_Complete/temp_${i}.fna.gz \| \
    sed \"s/^\>/\>${i}__/ \; s/ .*//\" \| \
    gzip \> Genomes_Complete/${i}.fna.gz \; \
    rm -f Genomes_Complete/temp_${i}.fna.gz
done > do_edit.sh
parallel -a do_edit.sh -j ${VAR_PARALLEL_EDIT}
rm -f do_edit.sh temp_filelist

find Genomes_Incomplete/ -type f -name '*.fna.gz' | sed 's/GTDB\/Genomes_Incomplete\///' | sed 's/\.fna\.gz//' > temp_filelist
for i in $(cat temp_filelist)
do
    echo mv Genomes_Incomplete/${i}.fna.gz Genomes_Incomplete/temp_${i}.fna.gz \; \
    zcat Genomes_Incomplete/temp_${i}.fna.gz \| \
    sed \"s/^\>/\>${i}__/ \; s/ .*//\" \| \
    gzip \> Genomes_Incomplete/${i}.fna.gz \; \
    rm -f Genomes_Incomplete/temp_${i}.fna.gz
done > do_edit.sh
parallel -a do_edit.sh -j ${VAR_PARALLEL}
rm -f do_edit.sh temp_filelist

# make directories to store downloaded annotation information
mkdir Annotation_Complete
mkdir Annotation_Incomplete

# download annotation information
awk 'BEGIN{FS=OFS="/";filesuffix="genomic.gff.gz"}{ftpdir=$0;asm=$10;file=asm"_"filesuffix;print "wget -P Annotation_Complete/ "ftpdir,file}' \
    ftpdirpaths_genbank_complete > ftpdownload
parallel -a ftpdownload -j ${VAR_PARALLEL_DL}
rm -f ftpdownload

awk 'BEGIN{FS=OFS="/";filesuffix="genomic.gff.gz"}{ftpdir=$0;asm=$10;file=asm"_"filesuffix;print "wget -P Annotation_Complete/ "ftpdir,file}' \
    ftpdirpaths_refseq_complete > ftpdownload
parallel -a ftpdownload -j ${VAR_PARALLEL_DL}
rm -f ftpdownload

awk 'BEGIN{FS=OFS="/";filesuffix="genomic.gff.gz"}{ftpdir=$0;asm=$10;file=asm"_"filesuffix;print "wget -P Annotation_Incomplete/ "ftpdir,file}' \
    ftpdirpaths_genbank_incomplete > ftpdownload
parallel -a ftpdownload -j ${VAR_PARALLEL_DL}
rm -f ftpdownload

awk 'BEGIN{FS=OFS="/";filesuffix="genomic.gff.gz"}{ftpdir=$0;asm=$10;file=asm"_"filesuffix;print "wget -P Annotation_Incomplete/ "ftpdir,file}' \
ftpdirpaths_refseq_incomplete > ftpdownload
parallel -a ftpdownload -j ${VAR_PARALLEL_DL}
rm -f ftpdownload

# simplify annotation filenames
find Annotation_Complete/ -type f -name '*.gff.gz' | \
    awk -F '.[0-9]_' '{print "mv "$0" "$1".gff.gz"}' > rename.sh
sh rename.sh
find Annotation_Incomplete/ -type f -name '*.gff.gz' | \
    awk -F '.[0-9]_' '{print "mv "$0" "$1".gff.gz"}' > rename.sh
sh rename.sh
rm -f rename.sh

# identify unannotated genomes
Rscript Scripts/genomes_missing_annotation.R

# identify rRNA genes in unannotated genomes
for i in $(cat missingannotation_complete.txt)
do
    zcat Genomes_Complete/${i}.fna.gz | \
    barrnap --threads ${VAR_BARRNAP_THREADS} \
        --reject 0.8 > Annotation_Complete/${i}_rrna.gff
done

for i in $(cat missingannotation_incomplete.txt)
do
    zcat Genomes_Incomplete/${i}.fna.gz | \
        barrnap --threads ${VAR_BARRNAP_THREADS} \
            --reject 0.8 > Annotation_Incomplete/${i}_rrna.gff
done

# extract rRNA genes and edit annotation seqids to contain RefSeq or GenBank accession and match fasta headers
for i in $(find Annotation_Complete/ -type f -name '*.gff.gz' | sed 's/GTDB\/Annotation_Complete\/// ; s/\.gff\.gz//')
do
    echo zcat Annotation_Complete/${i}.gff.gz \| \
    awk \'\$3==\"rRNA\"\' \| \
    sed \"s/^/${i}__/\" \> Annotation_Complete/${i}_rrna.gff
done > getRNA.sh
parallel -a getRNA.sh -j ${VAR_PARALLEL_EDIT}
rm -f getRNA.sh

for i in $(find Annotation_Incomplete/ -type f -name '*.gff.gz' | sed 's/GTDB\/Annotation_Incomplete\/// ; s/\.gff\.gz//')
do
    echo zcat Annotation_Incomplete/${i}.gff.gz \| \
    awk \'\$3==\"rRNA\"\' \| \
    sed \"s/^/${i}__/\" \> \
    Annotation_Incomplete/${i}_rrna.gff
done > getRNA.sh
parallel -a getRNA.sh -j ${VAR_PARALLEL_EDIT}
rm -f getRNA.sh

# make directories to store outputs
mkdir Outputs_Complete/
mkdir Outputs_Incomplete/

# make single multifasta files for all assemblies and annotations
find Genomes_Complete/ -type f -name "*.fna.gz" \
    -exec zcat {} + > Outputs_Complete/combined.fna
find Genomes_Incomplete/ -type f -name "*.fna.gz" \
    -exec zcat {} + > Outputs_Incomplete/combined.fna
find Annotation_Complete/ -type f -name "*_rrna.gff" \
    -exec cat {} + > Outputs_Complete/combined_rrna.gff
find Annotation_Incomplete/ -type f -name "*_rrna.gff" \
    -exec cat {} + > Outputs_Incomplete/combined_rrna.gff

# calculate the length of all sequences in mulitfasta
# this info will be used later to discard operons which cross the start/end boundary of the chromosome
awk '/^>/ {if (seqlen) print seqlen;print;seqlen=0;next} {seqlen+=length($0)}END{print seqlen}' \
    Outputs_Complete/combined.fna | \
    sed 's/\.[0-9] .*//' | \
    sed 's/^>//' | \
    paste - - > Outputs_Complete/seq_length.tsv

awk '/^>/ {if (seqlen) print seqlen;print;seqlen=0;next} {seqlen+=length($0)}END{print seqlen}' \
    Outputs_Incomplete/combined.fna | \
    sed 's/\.[0-9] .*//' | \
    sed 's/^>//' | \
    paste - - > Outputs_Incomplete/seq_length.tsv

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
    -fi Outputs_Complete/combined.fna \
    -bed Outputs_Complete/full_operons.gff \
    -fo Outputs_Complete/full_operons.fna \
    -s
rm -f Outputs_Complete/combined.fna.fai

bedtools getfasta \
    -fi Outputs_Incomplete/combined.fna \
    -bed Outputs_Incomplete/full_operons.gff \
    -fo Outputs_Incomplete/full_operons.fna \
    -s
rm -f Outputs_Incomplete/combined.fna.fai

# remove combined assembly and annotation files to save space
# temporarily keeping individual assembly files to do pairwise ANI
rm -f Outputs_Complete/combined.fna Outputs_Incomplete/combined.fna
rm -f Outputs_Complete/combined_rrna.gff Outputs_Combined/combined_rrna.gff

# renames sequences in operon fasta file by giving each operon a unique identifier
cut -f 2 Outputs_Complete/operon_identifiers.txt | \
    sed 's/^/>/' > Outputs_Complete/operon_identifiers_headers.txt
cut -f 2 Outputs_Incomplete/operon_identifiers.txt | \
    sed 's/^/>/' > Outputs_Incomplete/operon_identifiers_headers.txt

# replaces the fasta headers with unique operon identifiers
sed -i 's/>.* />/' Outputs_Complete/operon_identifiers_headers.txt
awk 'NR%2==0' Outputs_Complete/full_operons.fna | \
    paste -d '\n' Outputs_Complete/operon_identifiers_headers.txt - > Outputs_Complete/full_operons_identifiers.fna
sed -i 's/>.* />/' Outputs_Incomplete/operon_identifiers_headers.txt
awk 'NR%2==0' Outputs_Incomplete/full_operons.fna | \
    paste -d '\n' Outputs_Incomplete/operon_identifiers_headers.txt - > RefSeq/Outputs_Incomplete/full_operons_identifiers.fna

# sorts sequences by length (important for NR clustering, will retain longest sequence in each cluster)
sortbyname.sh \
    in=Outputs_Complete/full_operons_identifiers.fna \
    out=Outputs_Complete/full_operons_sorted.fna \
    length \
    descending
rm -f Outputs_Complete/full_operons_identifiers.fna

sortbyname.sh \
    in=Outputs_Incomplete/full_operons_identifiers.fna \
    out=Outputs_Incomplete/full_operons_sorted.fna \
    length \
    descending
rm -f Outputs_Incomplete/full_operons_identifiers.fna

# removes redundant sequences from Complete Assemblies dataset based on specified identity cutoff
# also extracts relevant cluster information
vsearch \
    --cluster_smallmem Outputs_Complete/full_operons_sorted.fna \
    --id ${VAR_CLUSTERID_VSEARCH} \
    --centroids Outputs_Complete/nrRep.fna \
    --uc Outputs_Complete/clusters.tsv \
    --threads ${VAR_THREADS_VSEARCH} \
    -log Outputs_Complete/vsearchlog.txt \
    --consout Outputs_Complete/nrCon.fna
cat Outputs_Complete/clusters.tsv | \
    awk '$1=="S"' > Outputs_Complete/vsearch_centroids.tsv
cat Outputs_Complete/clusters.tsv | \
    awk '$1=="C"' > Outputs_Complete/vsearch_clusters.tsv
cat Outputs_Complete/clusters.tsv | \
    awk '$1=="H"' > Outputs_Complete/vsearch_hits.tsv

# adds operon sequences from Incomplete Assemblies to the NR dataset generated from complete assemblies, reruns clustering, and extracts relevant cluster information
# only adds NR sequences from Incomplete genome assemblies if they are not already represented in the Complete dataset
# also combines the taxonomy tables of these incomplete operons
mkdir Outputs_Combined/
cat Outputs_Complete/nrRep.fna Outputs_Incomplete/full_operons_sorted.fna > Outputs_Combined/full_operons_sorted.fna
vsearch \
    --cluster_smallmem Outputs_Combined/full_operons_sorted.fna \
    --id ${VAR_CLUSTERID_VSEARCH} \
    --centroids Outputs_Combined/nrRep.fna \
    --uc Outputs_Combined/clusters.tsv \
    --threads ${VAR_THREADS_VSEARCH} \
    -log Outputs_Combined/vsearchlog.txt \
    --consout Outputs_Combined/nrCon.fna \
    --usersort
cat Outputs_Combined/clusters.tsv | \
    awk '$1=="S"' > Outputs_Combined/vsearch_centroids.tsv
cat Outputs_Combined/clusters.tsv | \
    awk '$1=="C"' > Outputs_Combined/vsearch_clusters.tsv
cat Outputs_Combined/clusters.tsv | \
    awk '$1=="H"' > Outputs_Combined/vsearch_hits.tsv

# edit fasta headers of nrCon sequence files
sed -i 's/centroid=// ; s/;seqs=.*//' Outputs_Complete/nrCon.fna 
sed -i 's/centroid=// ; s/;seqs=.*//' Outputs_Combined/nrCon.fna

# move database files to single directory
mkdir -p Database_GTDB/Full/
mkdir -p Database_GTDB/NR/

sortbyname.sh \
    in=Outputs_Combined/full_operons_sorted.fna \
    out=Database_GTDB/Full/full.fna.gz

sortbyname.sh \
    in=Outputs_Combined/nrRep.fna \
    out=Database_GTDB/NR/nrRep.fna.gz

sortbyname.sh \
    in=Outputs_Combined/nrCon.fna \
    out=Database_GTDB/NR/nrCon.fna.gz

# remove database sequence files from output directories
rm -f Outputs_Combined/full_operons_sorted.fna
rm -f Outputs_Combined/nrRep.fna
rm -f Outputs_Combined/nrCon.fna

# get GTDB taxonomy for each rrn operon in databases
Rscript Scripts/gtdb_get_taxonomy.R
mv Outputs_Combined/taxFull.tsv Database_GTDB/Full/taxFull.tsv
mv Outputs_Combined/taxRep.tsv Database_GTDB/NR/taxRep.tsv
mv Outputs_Combined/taxLCA.tsv Database_GTDB/NR/taxLCA.tsv
mv Outputs_Combined/taxMaj.tsv Database_GTDB/NR/taxMaj.tsv

