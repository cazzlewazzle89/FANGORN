#load Fangorn conda environment
source /home/cwwalsh/miniconda3/etc/profile.d/conda.sh
conda activate fastani

ls Fangorn_GTDB/Genomes_Complete/*.fna.gz > genomelist.txt
split -l 2000 genomelist.txt

for i in $(ls xa*)
do
    fastANI --ql ${i} --rl genomelist.txt -t 20 -o Fangorn_GTDB/fastani_complete_${i}.tsv
done

cat Fangorn_GTDB/fastani_complete_xa*.tsv > Fangorn_GTDB/fastani_complete.tsv
rm -f Fangorn_GTDB/fastani_complete_xa*.tsv
rm -f xa*

sed -i 's|Fangorn_GTDB\/Genomes_Complete\/||g ; s|\.fna\.gz||g' Fangorn_GTDB/fastani_complete.tsv

# deactivate conda environment
conda deactivate

mail -s 'FastANI Complete' calum.walsh@unimelb.edu.au < /dev/null
