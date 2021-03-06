# testing full-length primer sequences to see which pair gives the highest % of amplicons

for i in "RefSeq" "GTDB"
do
    for j in "Complete" "Combined"
    do
        for k in "full" "nrRep"
        do
        
        in_silico_PCR.pl \
            -s Database/${i}/${j}/${k}.fna.gz \
            -p primerpairs.txt \
            -l 10000 \
            > FilesForPlotting/${i}_primertest_description_${j}${k}.txt \
            2> FilesForPlotting/${i}_primertest_amplicons_${j}${k}.fna
        
        done
    done
done

rename 'RefSeq' 'refseq' FilesForPlotting/RefSeq_*
rename 'GTDB' 'gtdb' FilesForPlotting/*

# send email notification that the script is complete
mail -s 'InSilicoPCR Complete' calum.walsh@unimelb.edu.au < /dev/null
