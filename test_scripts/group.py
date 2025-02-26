import re

def parse_chromosomes(content):
    chromosomes = []
    for line in content.split('\n'):
        line = line.strip()
        if not line:
            continue
        match = re.match(r"Chromosome (.+): (\d+) reads", line)
        if match:
            name = f"Chromosome {match.group(1)}"
            reads = int(match.group(2))
            chromosomes.append((name, reads))
    return chromosomes

def group_chromosomes(chroms, k):
    groups = [{'total': 0, 'members': []} for _ in range(k)]
    for name, reads in chroms:
        min_group = min(groups, key=lambda g: g['total'])
        min_group['total'] += reads
        min_group['members'].append((name, reads))
    return groups

def main():
    content = """Chromosome 1: 5484012 reads
Chromosome 2: 3197500 reads
Chromosome 3: 2846064 reads
Chromosome 4: 1516919 reads
Chromosome 5: 3327662 reads
Chromosome 6: 2403035 reads
Chromosome 7: 2296350 reads
Chromosome 8: 1601507 reads
Chromosome 9: 2528177 reads
Chromosome 10: 1543144 reads
Chromosome 11: 2864396 reads
Chromosome 12: 3160297 reads
Chromosome 13: 713680 reads
Chromosome 14: 1439826 reads
Chromosome 15: 1509406 reads
Chromosome 16: 2276522 reads
Chromosome 17: 3315912 reads
Chromosome 18: 490533 reads
Chromosome 19: 2968499 reads
Chromosome 20: 1951910 reads
Chromosome 21: 844783 reads
Chromosome 22: 1092565 reads
Chromosome OTHER: 383759 reads
Chromosome MT: 789251 reads
Chromosome X: 1966954 reads
Chromosome Y: 59630 reads
"""
    
    chromosomes = parse_chromosomes(content)
    chromosomes.sort(key=lambda x: -x[1])
    
    for k in range(1, 57):
        groups = group_chromosomes(chromosomes, k)
        print(f"\n分为 {k} 组的分组情况（总reads数均衡）：")
        for i, group in enumerate(groups, 1):
            print(f"组 {i}: 总reads数 = {group['total']}")
            print(f"包含 {len(group['members'])} 条染色体")
            for name, reads in group['members']:
                print(f"  {name}: {reads} reads")

if __name__ == "__main__":
    main()
