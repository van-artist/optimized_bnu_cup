#!/bin/bash
set -eo pipefail
trap 'echo "错误发生在第 $LINENO 行，退出状态码 $?" >&2' ERR

# 环境变量配置
java_bin="/es01/paratera/sce3640/m5c/packages/java/jdk-17.0.12/bin/java"
htsjdk_jar="/es01/paratera/sce3640/m5c/packages/java/htsjdk-4.1.3-9-gc31bc92-SNAPSHOT.jar"
umicollapse_jar="/es01/paratera/sce3640/m5c/packages/UMICollapse-1.0.0/umicollapse.jar"
snappy_jar="/es01/paratera/sce3640/m5c/packages/java/snappy-java-1.1.9.1.jar"

hisat_3n_build="/es01/paratera/sce3640/m5c/packages/hisat3n/hisat-3n/hisat-3n-build"
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
OUTPUT_DIR="/es01/paratera/sce3640/m5c/final/single/workspace/origin"

REF_GENOME="/es01/paratera/sce3640/m5c/m5C-UBSseq-0.1/reference/Homo_sapiens.GRCh38.dna.primary_assembly.fa"
REF_NCRNA="/es01/paratera/sce3640/m5c/m5C-UBSseq-0.1/reference/Homo_sapiens.GRCh38.ncrna.fa"

SAMPLES=("SRR23538290" "SRR23538291" "SRR23538292")

mkdir -p $OUTPUT_DIR

start_time_total=$(date +%s)

# 定义样本处理函数
process_sample() {
    local SAMPLE_ID=$1
    
    echo "===== Processing $SAMPLE_ID ====="
    
    # 创建样本输出目录
    local SAMPLE_DIR="${OUTPUT_DIR}/${SAMPLE_ID}"
    mkdir -p $SAMPLE_DIR

    # STAGE 1处理流程
    # 1. 前处理
    echo "[${SAMPLE_ID}] 步骤1: 前处理"
    step_start=$(date +%s)
    $cutseq ${DATA_DIR}/${SAMPLE_ID}.fastq -t 56 -A INLINE -m 20 --trim-polyA --ensure-inline-barcode \
        -o ${SAMPLE_DIR}/${SAMPLE_ID}.fastq_cut \
        -s ${SAMPLE_DIR}/${SAMPLE_ID}.fastq_tooshort \
        -u ${SAMPLE_DIR}/${SAMPLE_ID}.fastq_untrimmed
    step_end=$(date +%s)
    step_duration=$((step_end - step_start))
    echo "[${SAMPLE_ID}] 步骤1完成，耗时 ${step_duration}秒"
    
    # 2. ncRNA比对
    echo "[${SAMPLE_ID}] 步骤2: ncRNA比对"
    step_start=$(date +%s)
    $hisat_3n --index $REF_NCRNA \
        --summary-file ${SAMPLE_DIR}/map2ncrna.output.summary \
        --new-summary -q -U ${SAMPLE_DIR}/${SAMPLE_ID}.fastq_cut \
        -p 56 --base-change C,T --mp 8,2 --no-spliced-alignment --directional-mapping | \
        $samtools view -@ 16 -e '!flag.unmap' -O BAM \
        -U ${SAMPLE_DIR}/${SAMPLE_ID}.ncrna.unmapped.bam \
        -o ${SAMPLE_DIR}/${SAMPLE_ID}.ncrna.mapped.bam
    step_end=$(date +%s)
    step_duration=$((step_end - step_start))
    echo "[${SAMPLE_ID}] 步骤2完成，耗时 ${step_duration}秒"
    
    # 3. 提取未比对reads
    echo "[${SAMPLE_ID}] 步骤3: 提取未比对reads"
    step_start=$(date +%s)
    $samtools fastq -@ 56 -O ${SAMPLE_DIR}/${SAMPLE_ID}.ncrna.unmapped.bam > ${SAMPLE_DIR}/${SAMPLE_ID}.mRNA.fastq
    step_end=$(date +%s)
    step_duration=$((step_end - step_start))
    echo "[${SAMPLE_ID}] 步骤3完成，耗时 ${step_duration}秒"

    # 4. 基因组比对
    echo "[${SAMPLE_ID}] 步骤4: 基因组比对"
    step_start=$(date +%s)
    $hisat_3n --index $REF_GENOME -p 56 \
        --summary-file ${SAMPLE_DIR}/map2genome.output.summary \
        --new-summary -q -U ${SAMPLE_DIR}/${SAMPLE_ID}.mRNA.fastq \
        --directional-mapping --base-change C,T --pen-noncansplice 20 --mp 4,1 | \
        $samtools view -@ 16 -e '!flag.unmap' -O BAM \
        -U ${SAMPLE_DIR}/${SAMPLE_ID}.mRNA.genome.unmapped.bam \
        -o ${SAMPLE_DIR}/${SAMPLE_ID}.mRNA.genome.mapped.bam
    step_end=$(date +%s)
    step_duration=$((step_end - step_start))
    echo "[${SAMPLE_ID}] 步骤4完成，耗时 ${step_duration}秒"

    # 5. 排序和统计
    echo "[${SAMPLE_ID}] 步骤5: 排序和统计"
    step_start=$(date +%s)
    $samtools sort -@ 56 --write-index -O BAM \
        -o ${SAMPLE_DIR}/${SAMPLE_ID}.mRNA.genome.mapped.sorted.bam \
        ${SAMPLE_DIR}/${SAMPLE_ID}.mRNA.genome.mapped.bam

    $samtools view -@ 56 -F 3980 -c ${SAMPLE_DIR}/${SAMPLE_ID}.mRNA.genome.mapped.sorted.bam > \
        ${SAMPLE_DIR}/${SAMPLE_ID}.mRNA.genome.mapped.sorted.bam.tsv
    step_end=$(date +%s)
    step_duration=$((step_end - step_start))
    echo "[${SAMPLE_ID}] 步骤5完成，耗时 ${step_duration}秒"


    # 6. UMI去重
    echo "[${SAMPLE_ID}] 步骤6: UMI去重"
    step_start=$(date +%s)
    $java_bin -server -Xms8G -Xmx40G -Xss100M -Djava.io.tmpdir=$SAMPLE_DIR \
        -cp "${htsjdk_jar}:${umicollapse_jar}:${snappy_jar}" umicollapse.main.Main bam \
        -t 2 -T 16 --data naive --merge avgqual --two-pass \
        -i ${SAMPLE_DIR}/${SAMPLE_ID}.mRNA.genome.mapped.sorted.bam \
        -o ${SAMPLE_DIR}/${SAMPLE_ID}.mRNA.genome.mapped.sorted.dedup.bam \
        > ${SAMPLE_DIR}/${SAMPLE_ID}.mRNA.genome.mapped.sorted.dedup.log

    step_end=$(date +%s)
    step_duration=$((step_end - step_start))
    echo "[${SAMPLE_ID}] 步骤6完成，耗时 ${step_duration}秒"
    
    # 7. 生成索引
    echo "[${SAMPLE_ID}] 步骤7: 生成索引"
    step_start=$(date +%s)
    $samtools index -@ 56 \
        ${SAMPLE_DIR}/${SAMPLE_ID}.mRNA.genome.mapped.sorted.dedup.bam \
        ${SAMPLE_DIR}/${SAMPLE_ID}.mRNA.genome.mapped.sorted.dedup.bam.bai
    step_end=$(date +%s)
    step_duration=$((step_end - step_start))
    echo "[${SAMPLE_ID}] 步骤7完成，耗时 ${step_duration}秒"
    
    # 8. 过滤BAM
    echo "[${SAMPLE_ID}] 步骤8: 过滤BAM"
    step_start=$(date +%s)
    $samtools view -@ 56 -e "[XM] * 20 <= (qlen-sclen) && [Zf] <= 3 && 3 * [Zf] <= [Zf] + [Yf]" \
        ${SAMPLE_DIR}/${SAMPLE_ID}.mRNA.genome.mapped.sorted.dedup.bam \
        -O BAM -o ${SAMPLE_DIR}/${SAMPLE_ID}.mRNA.genome.mapped.sorted.dedup.filtered.bam
    step_end=$(date +%s)
    step_duration=$((step_end - step_start))
    echo "[${SAMPLE_ID}] 步骤8完成，耗时 ${step_duration}秒"

    # 9. 生成TSV和过滤后TSV
    echo "[${SAMPLE_ID}] 步骤9: 生成TSV和过滤后TSV"
    step_start=$(date +%s)
    {
    $samtools view -e "rlen<100000" -h ${SAMPLE_DIR}/${SAMPLE_ID}.mRNA.genome.mapped.sorted.dedup.bam | \
        $hisat_3n_table -p 56 -u --alignments - --ref $REF_GENOME --output-name /dev/stdout --base-change C,T | \
        cut -f 1,2,3,5,7 | bgzip -c > ${SAMPLE_DIR}/${SAMPLE_ID}_unfiltered_uniq.tsv.gz

    $samtools view -e "rlen<100000" -h ${SAMPLE_DIR}/${SAMPLE_ID}.mRNA.genome.mapped.sorted.dedup.bam | \
        $hisat_3n_table -p 56 -m --alignments - --ref $REF_GENOME --output-name /dev/stdout --base-change C,T | \
        cut -f 1,2,3,5,7 | bgzip -c > ${SAMPLE_DIR}/${SAMPLE_ID}_unfiltered_multi.tsv.gz
    
    $samtools view -e "rlen<100000" -h ${SAMPLE_DIR}/${SAMPLE_ID}.mRNA.genome.mapped.sorted.dedup.filtered.bam | \
        $hisat_3n_table -p 56 -u --alignments - --ref $REF_GENOME --output-name /dev/stdout --base-change C,T | \
        cut -f 1,2,3,5,7 | bgzip -c > ${SAMPLE_DIR}/${SAMPLE_ID}_filtered_uniq.tsv.gz

    $samtools view -e "rlen<100000" -h ${SAMPLE_DIR}/${SAMPLE_ID}.mRNA.genome.mapped.sorted.dedup.filtered.bam | \
        $hisat_3n_table -p 56 -m --alignments - --ref $REF_GENOME --output-name /dev/stdout --base-change C,T | \
        cut -f 1,2,3,5,7 | bgzip -c > ${SAMPLE_DIR}/${SAMPLE_ID}_filtered_multi.tsv.gz
    }
    step_end=$(date +%s)
    step_duration=$((step_end - step_start))
    echo "[${SAMPLE_ID}] 步骤9完成，耗时 ${step_duration}秒"

    # 10. 合并结果
    echo "[${SAMPLE_ID}] 步骤10: 合并结果"
    step_start=$(date +%s)
    $join_pileup_py -i \
        ${SAMPLE_DIR}/${SAMPLE_ID}_unfiltered_uniq.tsv.gz \
        ${SAMPLE_DIR}/${SAMPLE_ID}_unfiltered_multi.tsv.gz \
        ${SAMPLE_DIR}/${SAMPLE_ID}_filtered_uniq.tsv.gz \
        ${SAMPLE_DIR}/${SAMPLE_ID}_filtered_multi.tsv.gz \
        -o ${SAMPLE_DIR}/${SAMPLE_ID}_genome.arrow
    step_end=$(date +%s)
    step_duration=$((step_end - step_start))
    echo "[${SAMPLE_ID}] 步骤10完成，耗时 ${step_duration}秒" 
}

# Stage 1：处理样本
echo "===== Stage 1: Sample Processing ====="
for SAMPLE in "${SAMPLES[@]}"; do
    process_sample $SAMPLE
done

# Stage 2：结果聚合
echo "===== Stage 2: Result Aggregation ====="

echo "步骤1: group_pileup"
step_start=$(date +%s)
$group_pileup_py -i \
    ${OUTPUT_DIR}/SRR23538290/SRR23538290_genome.arrow \
    ${OUTPUT_DIR}/SRR23538291/SRR23538291_genome.arrow \
    ${OUTPUT_DIR}/SRR23538292/SRR23538292_genome.arrow \
    -o ${OUTPUT_DIR}/WT.arrow
step_end=$(date +%s)
step_duration=$((step_end - step_start))
echo "步骤1完成，耗时 ${step_duration}秒" 


echo "步骤2: select_sites"
step_start=$(date +%s)
$select_sites_py -i ${OUTPUT_DIR}/WT.arrow -o ${OUTPUT_DIR}/WT.prefilter.tsv
step_end=$(date +%s)
step_duration=$((step_end - step_start))
echo "步骤2完成，耗时 ${step_duration}秒" 

# 各样本最终过滤
echo "步骤3: filter_sites"
step_start=$(date +%s)
for SAMPLE in "${SAMPLES[@]}"; do
    $filter_sites_py -i ${OUTPUT_DIR}/${SAMPLE}/${SAMPLE}_genome.arrow \
        -m ${OUTPUT_DIR}/WT.prefilter.tsv \
        -b ${OUTPUT_DIR}/${SAMPLE}/${SAMPLE}.bg.tsv \
        -o ${OUTPUT_DIR}/${SAMPLE}/${SAMPLE}.filtered.tsv
done
step_end=$(date +%s)
step_duration=$((step_end - step_start))
echo "步骤3完成，耗时 ${step_duration}秒" 

echo "======= All Processing Completed ======="

end_time_total=$(date +%s)
total_time=$((end_time_total - start_time_total))
echo "======= 全部处理完成，总运行时间: ${total_time} 秒 ======="