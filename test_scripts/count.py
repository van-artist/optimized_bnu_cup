def count_reads_per_chromosome(sam_file_path):
    chromosome_counts = {}
    
    with open(sam_file_path, 'r') as file:
        for line in file:
            if line.startswith('@'):
                continue
            fields = line.strip().split('\t')
            if len(fields) < 3:
                continue
            
            flag = int(fields[1])
            if (flag & 0x900) != 0:
                continue
            
            chromosome = fields[2].replace("chr", "")
            if chromosome != '*':
                chromosome_counts[chromosome] = chromosome_counts.get(chromosome, 0) + 1
    
    return chromosome_counts

if __name__ == "__main__":
    sam_file = "/es01/paratera/sce3640/m5c/m5C-UBSseq-0.1/shellspace/SRR23538290/SRR23538290_filtered.sam"
    counts = count_reads_per_chromosome(sam_file)
    
    for chromosome in sorted(counts.keys(), key=lambda x: (int(x) if x.isdigit() else float('inf'), x)):
        print(f"Chromosome {chromosome}: {counts[chromosome]} reads")
