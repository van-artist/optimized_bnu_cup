#!/bin/bash

# Data directories and package directories
RESULT_DIR="./result"
DATA_DIR="./data" 
REF_DIR="$DATA_DIR/ref"
RAW_DATA_DIR="$DATA_DIR/raw"
PACKAGES_DIR="./packages"
SCRIPTS_DIR="./scripts"
WORKSPACE_DIR="./workspace"
OUTPUT_DIR="$WORKSPACE_DIR"
JAVA_DIR="$PACKAGES_DIR/java"
HISAT3N_DIR="$PACKAGES_DIR/hisat3n"
M5C_UBSSEQ_DIR="$PACKAGES_DIR/m5C-UBSseq"
mkdir -p "$RESULT_DIR" "$WORKSPACE_DIR"

# HISAT3N compilation tool path
HISAT3N_BUILD="$PACKAGES_DIR/hisat-3n/hisat-3n-build"

GENOME_FA="$REF_DIR/Homo_sapiens.GRCh38.dna.primary_assembly.fa"
NCRNA_FA="$REF_DIR/Homo_sapiens.GRCh38.ncrna.fa"

# Index output paths
REF_GENOME="$REF_DIR/Homo_sapiens.GRCh38.dna.primary_assembly.fa"
REF_NCRNA="$REF_DIR/Homo_sapiens.GRCh38.ncrna.fa"

# Call external script with path parameters
bash "$SCRIPTS_DIR/build_ref.sh" "$HISAT3N_BUILD" "$GENOME_FA" "$REF_GENOME" "$NCRNA_FA" "$REF_NCRNA" >build_stdout.txt

JAVA_BASE_DIR="$PACKAGES_DIR/java"
java_bin="$JAVA_BASE_DIR/jdk-17.0.12/bin/java"
htsjdk_jar="$JAVA_BASE_DIR/htsjdk-4.1.3-9-gc31bc92-SNAPSHOT.jar"
umicollapse_jar="$JAVA_BASE_DIR//umicollapse.jar"
snappy_jar="$JAVA_BASE_DIR/snappy-java-1.1.9.1.jar"
chmod +x "$java_bin"

HISAT_BASE_DIR="$PACKAGES_DIR/hisat3n"
hisat_3n_build="$HISAT_BASE_DIR/hisat-3n-build"
hisat_3n="$HISAT_BASE_DIR/hisat-3n"
hisat_3n_table="$HISAT_BASE_DIR/hisat-3n-table"

cutseq="cutseq"
samtools="samtools"
bgzip="bgzip"

PY_BIN_BASE_DIR="$M5C_UBSSEQ_DIR/bin"
join_pileup_py="$PY_BIN_BASE_DIR/join_pileup.py"
group_pileup_py="$PY_BIN_BASE_DIR/group_pileup.py"
select_sites_py="$PY_BIN_BASE_DIR/select_sites.py"
filter_sites_py="$PY_BIN_BASE_DIR/filter_sites.py"

chmod +x "$PY_BIN_BASE_DIR/join_pileup.py" \
         "$PY_BIN_BASE_DIR/group_pileup.py" \
         "$PY_BIN_BASE_DIR/select_sites.py" \
         "$PY_BIN_BASE_DIR/filter_sites.py"


SAMPLES=("SRR23538290" "SRR23538291" "SRR23538292")

start_time_total=$(date +%s)

mkdir -p $OUTPUT_DIR
start_time_total=$(date +%s)

# Define sample processing function
process_sample() {
    local SAMPLE_ID=$1
    
    echo "===== Processing $SAMPLE_ID ====="
    
    # Create sample output directory
    local SAMPLE_DIR="${OUTPUT_DIR}/${SAMPLE_ID}"
    mkdir -p $SAMPLE_DIR

    # STAGE 1 processing pipeline
    # 1. Preprocessing
    echo "[${SAMPLE_ID}] Step 1: Preprocessing"
    step_start=$(date +%s)
    $cutseq ${RAW_DATA_DIR}/${SAMPLE_ID}.fastq -t 20 -A INLINE -m 20 --trim-polyA --ensure-inline-barcode \
        -o ${SAMPLE_DIR}/${SAMPLE_ID}.fastq_cut \
        -s ${SAMPLE_DIR}/${SAMPLE_ID}.fastq_tooshort \
        -u ${SAMPLE_DIR}/${SAMPLE_ID}.fastq_untrimmedxxx
    step_end=$(date +%s)
    step_duration=$((step_end - step_start))
    echo "[${SAMPLE_ID}] Step 1 completed in ${step_duration} seconds"
    
    # 2. ncRNA alignment
    echo "[${SAMPLE_ID}] Step 2: ncRNA alignment"
    step_start=$(date +%s)
    $hisat_3n --index $REF_NCRNA \
        --summary-file ${SAMPLE_DIR}/map2ncrna.output.summary \
        --new-summary -q -U ${SAMPLE_DIR}/${SAMPLE_ID}.fastq_cut \
        -p 16 --base-change C,T --mp 8,2 --no-spliced-alignment --directional-mapping | \
        $samtools view -@ 16 -e '!flag.unmap' -O BAM \
        -U ${SAMPLE_DIR}/${SAMPLE_ID}.ncrna.unmapped.bam \
        -o ${SAMPLE_DIR}/${SAMPLE_ID}.ncrna.mapped.bam
    step_end=$(date +%s)
    step_duration=$((step_end - step_start))
    echo "[${SAMPLE_ID}] Step 2 completed in ${step_duration} seconds"
    
    # 3. Extract unmapped reads
    echo "[${SAMPLE_ID}] Step 3: Extract unmapped reads"
    step_start=$(date +%s)
    $samtools fastq -@ 16 -O ${SAMPLE_DIR}/${SAMPLE_ID}.ncrna.unmapped.bam > ${SAMPLE_DIR}/${SAMPLE_ID}.mRNA.fastq
    step_end=$(date +%s)
    step_duration=$((step_end - step_start))
    echo "[${SAMPLE_ID}] Step 3 completed in ${step_duration} seconds"

    # 4. Genome alignment
    echo "[${SAMPLE_ID}] Step 4: Genome alignment"
    step_start=$(date +%s)
    $hisat_3n --index $REF_GENOME -p 16 \
        --summary-file ${SAMPLE_DIR}/map2genome.output.summary \
        --new-summary -q -U ${SAMPLE_DIR}/${SAMPLE_ID}.mRNA.fastq \
        --directional-mapping --base-change C,T --pen-noncansplice 20 --mp 4,1 | \
        $samtools view -@ 16 -e '!flag.unmap' -O BAM \
        -U ${SAMPLE_DIR}/${SAMPLE_ID}.mRNA.genome.unmapped.bam \
        -o ${SAMPLE_DIR}/${SAMPLE_ID}.mRNA.genome.mapped.bam
    step_end=$(date +%s)
    step_duration=$((step_end - step_start))
    echo "[${SAMPLE_ID}] Step 4 completed in ${step_duration} seconds"

    # 5. Sorting and statistics
    echo "[${SAMPLE_ID}] Step 5: Sorting and statistics"
    step_start=$(date +%s)
    $samtools sort -@ 16 --write-index -O BAM \
        -o ${SAMPLE_DIR}/${SAMPLE_ID}.mRNA.genome.mapped.sorted.bam \
        ${SAMPLE_DIR}/${SAMPLE_ID}.mRNA.genome.mapped.bam

    $samtools view -@ 20 -F 3980 -c ${SAMPLE_DIR}/${SAMPLE_ID}.mRNA.genome.mapped.sorted.bam > \
        ${SAMPLE_DIR}/${SAMPLE_ID}.mRNA.genome.mapped.sorted.bam.tsv

    step_end=$(date +%s)
    step_duration=$((step_end - step_start))
    echo "[${SAMPLE_ID}] Step 5 completed in ${step_duration} seconds"

    # 6. UMI deduplication
    echo "[${SAMPLE_ID}] Step 6: UMI deduplication"
    step_start=$(date +%s)
    $java_bin -server -Xms16G -Xmx60G -Xss100M -XX:+UseG1GC -Djava.io.tmpdir=/dev/shm \
        -cp "${htsjdk_jar}:${umicollapse_jar}:${snappy_jar}" umicollapse.main.Main bam \
        --algo dir --data bktree --merge avgqual --two-pass \
        -i ${SAMPLE_DIR}/${SAMPLE_ID}.mRNA.genome.mapped.sorted.bam \
        -o ${SAMPLE_DIR}/${SAMPLE_ID}.mRNA.genome.mapped.sorted.dedup.bam \
        > ${SAMPLE_DIR}/${SAMPLE_ID}.mRNA.genome.mapped.sorted.dedup.log
    step_end=$(date +%s)
    step_duration=$((step_end - step_start))
    echo "[${SAMPLE_ID}] Step 6 completed in ${step_duration} seconds"
    
    # 7. Generate index
    echo "[${SAMPLE_ID}] Step 7: Generate index"
    step_start=$(date +%s)
    $samtools index -@ 8 \
        ${SAMPLE_DIR}/${SAMPLE_ID}.mRNA.genome.mapped.sorted.dedup.bam \
        ${SAMPLE_DIR}/${SAMPLE_ID}.mRNA.genome.mapped.sorted.dedup.bam.bai
    step_end=$(date +%s)
    step_duration=$((step_end - step_start))
    echo "[${SAMPLE_ID}] Step 7 completed in ${step_duration} seconds"
    
    # 8. Filter BAM
    echo "[${SAMPLE_ID}] Step 8: Filter BAM"
    step_start=$(date +%s)
    $samtools view -@ 8 -e "[XM] * 20 <= (qlen-sclen) && [Zf] <= 3 && 3 * [Zf] <= [Zf] + [Yf]" \
        ${SAMPLE_DIR}/${SAMPLE_ID}.mRNA.genome.mapped.sorted.dedup.bam \
        -O BAM -o ${SAMPLE_DIR}/${SAMPLE_ID}.mRNA.genome.mapped.sorted.dedup.filtered.bam
    step_end=$(date +%s)
    step_duration=$((step_end - step_start))
    echo "[${SAMPLE_ID}] Step 8 completed in ${step_duration} seconds"

    # 9. Convert BAM to SAM and split by chromosome
    echo "[${SAMPLE_ID}] Step 9: Convert BAM to SAM and split by chromosome"
    step_start=$(date +%s)
    {
        mkdir -p "${SAMPLE_DIR}/split"
        rm -rf "${SAMPLE_DIR}/split"/*

        $samtools view -@ 8 -e "rlen<100000" "${SAMPLE_DIR}/${SAMPLE_ID}.mRNA.genome.mapped.sorted.dedup.bam" > "${SAMPLE_DIR}/${SAMPLE_ID}_unfiltered.sam" &
        $samtools view -@ 8 -e "rlen<100000" "${SAMPLE_DIR}/${SAMPLE_ID}.mRNA.genome.mapped.sorted.dedup.filtered.bam" > "${SAMPLE_DIR}/${SAMPLE_ID}_filtered.sam" &
        wait

        declare -A CHROM_GROUP_CONFIG=(
            [2]="5 17 12 11 9 7 X 8 10 15 21 MT 18 Y|_OTHER_"
            [3]="1 19 9 16 8 15 MT 18|17 2 3 7 X 10 14 13 Y|_OTHER_"
            [4]="5 11 9 8 10 22|17 20 6 14 21 Y|2 12 7 X 15 13 18|19 19 3 7 15|_OTHER_"
            [5]="1 6 8 22|5 9 7 15 13 Y|17 3 20 4 21|12 19 X 10 MT|_OTHER_"
            [6]="1 16 13 Y|5 7 X 21 18|2 9 8 15|12 3 4 22|19 11 10 14|_OTHER_"
            [7]="1 8 Y|5 X 20|17 16 14 13|2 7 10 18|19 9 4 MT|11 3 22 21|_OTHER_"
            [8]="1 13 18|5 8 10|17 20 22|2 X 4|12 16 MT Y|19 7 15|11 6 14|_OTHER_"
            [9]="1 Y|5 4 14|17 10 22|2 8 15|12 20 13|19 X MT|11 16 18|9 6 21|_OTHER_"
            [10]="1|5 14 18|17 15|2 4 13|12 10 MT|19 8 22|11 20|3 X Y|6 7 21|_OTHER_"
        )

        process_sam() {
            local split=7
            local suffix=$1
            local sam_fil="${SAMPLE_DIR}/${SAMPLE_ID}_${suffix}.sam"
            
            local config_str="${CHROM_GROUP_CONFIG[$split]}"
            IFS='|' read -ra groups <<< "$config_str"
            local num_groups=${#groups[@]}

            for ((group_index=0; group_index<num_groups-1; group_index++)); do
                echo "${groups[$group_index]}" | tr ' ' '\n' > "${SAMPLE_DIR}/split/group${group_index}_${suffix}.chroms"
            done

            local combined="${SAMPLE_DIR}/split/combined_${suffix}.chroms"
            cat "${SAMPLE_DIR}"/split/group[0-9]*_${suffix}.chroms > "$combined"

            for ((group_index=0; group_index<num_groups; group_index++)); do                
                if [[ $group_index -eq $((num_groups-1)) ]] && [[ "${groups[$group_index]}" == *"_OTHER_"* ]]; then
                    awk -v combined="$combined" '
                        BEGIN {
                            while (getline line < combined > 0) { chroms[line]++ }
                        }
                        !($3 in chroms) { print $0 }
                    ' "$sam_fil" > "${SAMPLE_DIR}/split/${SAMPLE_ID}_${suffix}_part_${group_index}.sam" &
                else
                    local chrom_file="${SAMPLE_DIR}/split/group${group_index}_${suffix}.chroms"
                    awk '
                        NR==FNR { chroms[$1]++; next }
                        $3 in chroms { print $0 }
                    ' "$chrom_file" "$sam_fil" > "${SAMPLE_DIR}/split/${SAMPLE_ID}_${suffix}_part_${group_index}.sam" &
                fi
            done
            wait

            rm -f "${SAMPLE_DIR}/split/combined_${suffix}.chroms"
            rm -f "${SAMPLE_DIR}/split/group"*"_${suffix}.chroms"
        }

        process_sam "unfiltered" &
        process_sam "filtered" &
        wait
    }
    step_end=$(date +%s)
    echo "[${SAMPLE_ID}] Step 9 completed in $((step_end - step_start)) seconds"

    # 10. Generate TSV files
    echo "[${SAMPLE_ID}] Step 10: Generate TSV files"
    step_start=$(date +%s)
    {
        generate_tsv() {
            local suffix=$1
            for file in "${SAMPLE_DIR}/split/${SAMPLE_ID}_${suffix}_part_"*".sam"; do
                base_name="${file%.sam}"
                $hisat_3n_table -p 1 -u --alignments "$file" --ref "$REF_GENOME" \
                    --output-name "${base_name}_uniq.tsv" --base-change C,T &
                $hisat_3n_table -p 1 -m --alignments "$file" --ref "$REF_GENOME" \
                    --output-name "${base_name}_multi.tsv" --base-change C,T &
            done
            wait
        }
        
        generate_tsv "unfiltered" &
        generate_tsv "filtered" &
        wait
    }
    step_end=$(date +%s)
    echo "[${SAMPLE_ID}] Step 10 completed in $((step_end - step_start)) seconds"

    # 11. Merge TSV files
    echo "[${SAMPLE_ID}] Step 11: Merge TSV files"
    step_start=$(date +%s)
    {
        merge_tsv() {
            local suffix=$1
            local output_uniq="${SAMPLE_DIR}/${SAMPLE_ID}_${suffix}_uniq.tsv"
            local output_multi="${SAMPLE_DIR}/${SAMPLE_ID}_${suffix}_multi.tsv"
            
            awk 'FNR==1 && NR!=1 { next } 1' "${SAMPLE_DIR}/split/${SAMPLE_ID}_${suffix}_part_"*"_uniq.tsv" > "${output_uniq}" &
            awk 'FNR==1 && NR!=1 { next } 1' "${SAMPLE_DIR}/split/${SAMPLE_ID}_${suffix}_part_"*"_multi.tsv" > "${output_multi}" &
            wait
        }

        merge_tsv "unfiltered" &
        merge_tsv "filtered" &
        wait
    }
    rm -rf "${SAMPLE_DIR}/split"s
    step_end=$(date +%s)
    echo "[${SAMPLE_ID}] Step 11 completed in $((step_end - step_start)) seconds"
    
    # 12. Compress results
    echo "[${SAMPLE_ID}] Step 12: Compress results"
    step_start=$(date +%s)
    {
        cut -f 1,2,3,5,7 ${SAMPLE_DIR}/${SAMPLE_ID}_unfiltered_uniq.tsv | bgzip -@ 4 -c > ${SAMPLE_DIR}/${SAMPLE_ID}_unfiltered_uniq.tsv.gz &
        cut -f 1,2,3,5,7 ${SAMPLE_DIR}/${SAMPLE_ID}_unfiltered_multi.tsv | bgzip -@ 4 -c > ${SAMPLE_DIR}/${SAMPLE_ID}_unfiltered_multi.tsv.gz &
        cut -f 1,2,3,5,7 ${SAMPLE_DIR}/${SAMPLE_ID}_filtered_uniq.tsv | bgzip -@ 4 -c > ${SAMPLE_DIR}/${SAMPLE_ID}_filtered_uniq.tsv.gz &
        cut -f 1,2,3,5,7 ${SAMPLE_DIR}/${SAMPLE_ID}_filtered_multi.tsv | bgzip -@ 4 -c > ${SAMPLE_DIR}/${SAMPLE_ID}_filtered_multi.tsv.gz &
        wait
    }
    step_end=$(date +%s)
    echo "[${SAMPLE_ID}] Step 12 completed in $((step_end - step_start)) seconds"

    # 13. Merge final results
    echo "[${SAMPLE_ID}] Step 13: Merge final results"
    step_start=$(date +%s)
    $join_pileup_py -i \
        ${SAMPLE_DIR}/${SAMPLE_ID}_unfiltered_uniq.tsv.gz \
        ${SAMPLE_DIR}/${SAMPLE_ID}_unfiltered_multi.tsv.gz \
        ${SAMPLE_DIR}/${SAMPLE_ID}_filtered_uniq.tsv.gz \
        ${SAMPLE_DIR}/${SAMPLE_ID}_filtered_multi.tsv.gz \
        -o ${SAMPLE_DIR}/${SAMPLE_ID}_genome.arrow
    step_end=$(date +%s)
    echo "[${SAMPLE_ID}] Step 13 completed in $((step_end - step_start)) seconds"
}

# Stage 1: Parallel sample processing
echo "===== Stage 1: Parallel Sample Processing ====="
for SAMPLE in "${SAMPLES[@]}"; do
    process_sample $SAMPLE &
done
wait

# Stage 2: Result aggregation
echo "===== Stage 2: Result Aggregation ====="

echo "Step 1: group_pileup"
step_start=$(date +%s)
$group_pileup_py -i \
    ${OUTPUT_DIR}/SRR23538290/SRR23538290_genome.arrow \
    ${OUTPUT_DIR}/SRR23538291/SRR23538291_genome.arrow \
    ${OUTPUT_DIR}/SRR23538292/SRR23538292_genome.arrow \
    -o ${OUTPUT_DIR}/WT.arrow
step_end=$(date +%s)
step_duration=$((step_end - step_start))
echo "Step 1 completed in ${step_duration} seconds" 

echo "Step 2: select_sites"
step_start=$(date +%s)
$select_sites_py -i ${OUTPUT_DIR}/WT.arrow -o ${OUTPUT_DIR}/WT.prefilter.tsv
step_end=$(date +%s)
step_duration=$((step_end - step_start))
echo "Step 2 completed in ${step_duration} seconds" 

# Final filtering for each sample
echo "Step 3: filter_sites"
step_start=$(date +%s)
for SAMPLE in "${SAMPLES[@]}"; do
    $filter_sites_py -i ${OUTPUT_DIR}/${SAMPLE}/${SAMPLE}_genome.arrow \
        -m ${OUTPUT_DIR}/WT.prefilter.tsv \
        -b ${OUTPUT_DIR}/${SAMPLE}/${SAMPLE}.bg.tsv \
        -o ${RESULT_DIR}/${SAMPLE}.filtered.tsv
done
step_end=$(date +%s)
step_duration=$((step_end - step_start))
echo "Step 3 completed in ${step_duration} seconds" 

echo "======= All Processing Completed ======="

end_time_total=$(date +%s)
total_time=$((end_time_total - start_time_total))
echo "======= Total runtime: ${total_time} seconds ======="