#!/bin/bash
process_sample() {
    {
    # Process unfiltered unique reads (example section shown in detail)
    samtools view -e "rlen<100000" -h ${SAMPLE_DIR}/${SAMPLE_ID}.mRNA... |
    hisat_3n_table -p 1 -u --alignments - --ref $REF_GENOME... |
    cut -f 1,2,3,5,7 | bgzip -c > ${SAMPLE_DIR}/${SAMPLE_ID}_unfiltered_uniq.tsv.gz &
    # [Similar processing for 3 other combinations:
    # unfiltered_multi (-m) 
    # filtered_uniq (-u) or filtered_multi (-m) + filtered.bam]
    wait
    }
}
for SAMPLE in "${SAMPLES[@]}"; do
    process_sample $SAMPLE &
done
wait
