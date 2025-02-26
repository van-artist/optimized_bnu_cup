#!/bin/bash
#SBATCH --job-name=umi_test
#SBATCH --output=./out/umi_test/%x_%j.out
#SBATCH --error=./out/umi_test/%x_%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=56
#SBATCH --partition=G1Part_sce
#SBATCH --nodes=1

java_bin="/es01/paratera/sce3640/m5c/packages/java/jdk-17.0.12/bin/java"
htsjdk_jar="/es01/paratera/sce3640/m5c/packages/java/htsjdk-4.1.3-9-gc31bc92-SNAPSHOT.jar"
umicollapse_jar="/es01/paratera/sce3640/m5c/packages/UMICollapse-1.0.0/umicollapse.jar"
snappy_jar="/es01/paratera/sce3640/m5c/packages/java/snappy-java-1.1.9.1.jar"

mkdir -p umi_test

run_umicollapse() {
    local algo_desc="$1"
    local data_type="$2"
    local two_pass="$3"
    
    local suffix="${data_type}"
    [ -n "$two_pass" ] && suffix+="_twopass"
    local output_bam="/es01/paratera/sce3640/m5c/m5C-UBSseq-0.1/umi_test/umi_test_${suffix}.bam"
    local log_file="/es01/paratera/sce3640/m5c/m5C-UBSseq-0.1/umi_test/umi_test_${suffix}.log"
    local java_log_pipe=$(mktemp -u)

    # 创建日志管道
    mkfifo "$java_log_pipe" || { echo "无法创建命名管道"; exit 1; }
    tee "$log_file" < "$java_log_pipe" &
    
    echo "UMI去重 ${algo_desc}"
    step_start=$(date +%s)
    
    # 启动Java进程
    LC_ALL=C $java_bin -server -Xms24G -Xmx160G -Xss100M -XX:+UseG1GC -Djava.io.tmpdir=/dev/shm \
        -cp "${htsjdk_jar}:${umicollapse_jar}:${snappy_jar}" umicollapse.main.Main bam \
        --algo dir --data $data_type --merge avgqual $two_pass \
        -i /es01/paratera/sce3640/m5c/m5C-UBSseq-0.1/baseline/SRR23538291/SRR23538291.mRNA.genome.mapped.sorted.bam \
        -o $output_bam \
        > "$java_log_pipe" 2>&1 &
    java_pid=$!

    # 内存监控循环
    max_rss_kb=0
    while kill -0 $java_pid 2>/dev/null; do
        current_rss=$(ps -o rss= -p $java_pid | awk '{print $1}')
        if [ -n "$current_rss" ] && [ $current_rss -gt $max_rss_kb ]; then
            max_rss_kb=$current_rss
        fi
        sleep 1
    done

    # 等待进程结束获取状态码
    wait $java_pid
    java_exit_code=$?

    step_end=$(date +%s)
    step_duration=$((step_end - step_start))
    
    # 处理内存数据
    if [ $max_rss_kb -gt 0 ]; then
        max_rss_gb=$(echo "scale=2; $max_rss_kb / 1024 / 1024" | bc)
        mem_info="峰值内存 ${max_rss_gb}GB"
    else
        mem_info="内存监控失败"
    fi
    
    # 输出结果
    echo "耗时 ${step_duration}秒 | ${mem_info} | Java状态码 ${java_exit_code} | 输出文件 ${output_bam}"
    rm -f "$java_log_pipe"
    echo "----------------------------------------"
}

# 执行所有算法变体
run_umicollapse "bktree + two-pass" bktree "--two-pass"
run_umicollapse "fenwickbktree + two-pass" fenwickbktree "--two-pass"
run_umicollapse "ngrambktree + two-pass" ngrambktree "--two-pass"
run_umicollapse "bktree" bktree ""
run_umicollapse "fenwickbktree" fenwickbktree ""
run_umicollapse "ngrambktree" ngrambktree ""
