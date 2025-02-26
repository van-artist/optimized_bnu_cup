#!/bin/bash
set -eo pipefail
trap 'echo "错误发生在第 $LINENO 行，退出状态码 $?" >&2' ERR

split=5
while [[ "$#" -gt 0 ]]; do
    case $1 in
        -split) split="$2"; shift ;;
        *) echo "未知参数: $1"; exit 1 ;;
    esac
    shift
done

# 环境变量配置
java_bin="/es01/paratera/sce3640/m5c/packages/java/jdk-17.0.12/bin/java"
htsjdk_jar="/es01/paratera/sce3640/m5c/packages/java/htsjdk-4.1.3-9-gc31bc92-SNAPSHOT.jar"
umicollapse_jar="/es01/paratera/sce3640/m5c/packages/UMICollapse-1.0.0/umicollapse.jar"
snappy_jar="/es01/paratera/sce3640/m5c/packages/java/snappy-java-1.1.9.1.jar"

hisat_3n_build="/es01/paratera/sce3640/m5c/packages/hisat-3n/hisat-3n-build"
hisat_3n="/es01/paratera/sce3640/m5c/packages/hisat3n/hisat-3n/hisat-3n"
hisat_3n_table="/es01/paratera/sce3640/m5c/packages/hisat3n/hisat-3n/hisat-3n-table"

cutseq="cutseq"
samtools="samtools"
bgzip="bgzip"

join_pileup_py="/es01/paratera/sce3640/m5c/m5C-UBSseq-0.1/bin/join_pileup.py"
group_pileup_py="/es01/paratera/sce3640/m5c/m5C-UBSseq-0.1/bin/group_pileup.py"
select_sites_py="/es01/paratera/sce3640/m5c/m5C-UBSseq-0.1/bin/select_sites.py"
filter_sites_py="/es01/paratera/sce3640/m5c/m5C-UBSseq-0.1/bin/filter_sites.py"

DATA_DIR="/es01/paratera/sce3640/m5c/m5C-UBSseq-0.1/raw_data"
OUTPUT_DIR="/es01/paratera/sce3640/m5c/m5C-UBSseq-0.1/shelltablesplit/shelltable${split}"
BASELINE_DIR="/es01/paratera/sce3640/m5c/m5C-UBSseq-0.1/baseline"

REF_GENOME="/es01/paratera/sce3640/m5c/m5C-UBSseq-0.1/reference/Homo_sapiens.GRCh38.dna.primary_assembly.fa"
REF_NCRNA="/es01/paratera/sce3640/m5c/m5C-UBSseq-0.1/reference/Homo_sapiens.GRCh38.ncrna.fa"

SAMPLES=("SRR23538290" "SRR23538291" "SRR23538292")

mkdir -p $OUTPUT_DIR

start_time_total=$(date +%s)

# 定义样本处理函数
process_sample() {
    local SAMPLE_ID=$1
    
    echo "===== 分${split}组 Processing $SAMPLE_ID ====="
    
    # 创建样本输出目录
    local SAMPLE_DIR="${OUTPUT_DIR}/${SAMPLE_ID}"
    local INPUT_DIR="${BASELINE_DIR}/${SAMPLE_ID}"
    mkdir -p $SAMPLE_DIR

    # STAGE 1处理流程
    # 7. 生成索引
    echo "分${split}组 [${SAMPLE_ID}] 步骤7: 生成索引"
    step_start=$(date +%s)
    $samtools index -@ 8 \
        ${INPUT_DIR}/${SAMPLE_ID}.mRNA.genome.mapped.sorted.dedup.bam \
        ${SAMPLE_DIR}/${SAMPLE_ID}.mRNA.genome.mapped.sorted.dedup.bam.bai
    step_end=$(date +%s)
    step_duration=$((step_end - step_start))
    echo "分${split}组 [${SAMPLE_ID}] 步骤7完成，耗时 ${step_duration}秒"

    # 9. BAM转SAM并按染色体组分割
    echo "分${split}组 [${SAMPLE_ID}] 步骤9: BAM转SAM并按染色体组分割"
    step_start=$(date +%s)
    {
        mkdir -p "${SAMPLE_DIR}/split"
        rm -rf "${SAMPLE_DIR}/split"/*

        $samtools view -@ 8 -e "rlen<100000" "${INPUT_DIR}/${SAMPLE_ID}.mRNA.genome.mapped.sorted.dedup.bam" > "${SAMPLE_DIR}/${SAMPLE_ID}_unfiltered.sam" &
        $samtools view -@ 8 -e "rlen<100000" "${INPUT_DIR}/${SAMPLE_ID}.mRNA.genome.mapped.sorted.dedup.filtered.bam" > "${SAMPLE_DIR}/${SAMPLE_ID}_filtered.sam" &
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

            # rm -f "${SAMPLE_DIR}/split/group"*"_${suffix}.chroms"
        }

        process_sam "unfiltered" &
        process_sam "filtered" &
        wait
    }
    step_end=$(date +%s)
    echo "分${split}组 [${SAMPLE_ID}] 步骤9完成，耗时 $((step_end - step_start))秒"

    # 10. 生成TSV文件
    echo "分${split}组 [${SAMPLE_ID}] 步骤10: 生成TSV文件"
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
    echo "分${split}组 [${SAMPLE_ID}] 步骤10完成，耗时 $((step_end - step_start))秒"

    # 11. 合并TSV文件
    echo "分${split}组 [${SAMPLE_ID}] 步骤11: 合并TSV文件"
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
    step_end=$(date +%s)
    echo "分${split}组 [${SAMPLE_ID}] 步骤11完成，耗时 $((step_end - step_start))秒"

    # 12. 结果压缩
    echo "分${split}组 [${SAMPLE_ID}] 步骤12: 结果压缩"
    step_start=$(date +%s)
    {
        cut -f 1,2,3,5,7 ${SAMPLE_DIR}/${SAMPLE_ID}_unfiltered_uniq.tsv | bgzip -@ 4 -c > ${SAMPLE_DIR}/${SAMPLE_ID}_unfiltered_uniq.tsv.gz &
        cut -f 1,2,3,5,7 ${SAMPLE_DIR}/${SAMPLE_ID}_unfiltered_multi.tsv | bgzip -@ 4 -c > ${SAMPLE_DIR}/${SAMPLE_ID}_unfiltered_multi.tsv.gz &
        cut -f 1,2,3,5,7 ${SAMPLE_DIR}/${SAMPLE_ID}_filtered_uniq.tsv | bgzip -@ 4 -c > ${SAMPLE_DIR}/${SAMPLE_ID}_filtered_uniq.tsv.gz &
        cut -f 1,2,3,5,7 ${SAMPLE_DIR}/${SAMPLE_ID}_filtered_multi.tsv | bgzip -@ 4 -c > ${SAMPLE_DIR}/${SAMPLE_ID}_filtered_multi.tsv.gz &
        wait
    }
    step_end=$(date +%s)
    echo "分${split}组 [${SAMPLE_ID}] 步骤12完成，耗时 $((step_end - step_start))秒"

}

# Stage 1：并行处理样本
echo "===== 分${split}组 Stage 1: Parallel Sample Processing ====="
for SAMPLE in "${SAMPLES[@]}"; do
    process_sample $SAMPLE &
done
wait

echo "======= 分${split}组 All Processing Completed ======="

end_time_total=$(date +%s)
total_time=$((end_time_total - start_time_total))
echo "======= 分${split}组 全部处理完成，总运行时间: ${total_time} 秒 ======="