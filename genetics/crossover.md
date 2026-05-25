# Crossover

染色体互换（chromosomal crossover）是指有性生殖中，两条同源染色体的非姊妹染色单体部分的 DNA 发生交换的过程，是基因重组的最后阶段之一，发生在减数分裂前期 I 的粗线期的联会的过程中。联会在联会复合体发育前就开始了，直到前期 I 接近结束时才完成。交叉互换通常发生在配对染色体上的匹配区域断裂，然后与另一个染色体重新连接时。

染色体互换的起源有两种流行且重叠的理论，这两个理论来自关于减数分裂起源的不同理论。
第一个理论基于减数分裂进化为另一种 DNA 修复方法的想法，认为交叉互换是一种替代可能受损的 DNA 部分的新方法。
第二个理论来自减数分裂从细菌转化而来的想法变换，认为交叉互换具有产生多样性的功能。

I 类交叉（Class I crossovers）是由 ZMM 蛋白家族（如酵母中的 Zip1-4、Msh4/5、Mer3 及其同源蛋白）介导的主要交叉类型，通过形成并特异性解离双霍利迪连接体（dHJ）完成重组。其最显著的特征是存在交叉干扰，即一个位点发生交叉会降低邻近位点再次交叉的概率，从而确保交叉事件在染色体上均匀分布，有助于减数分裂过程中染色体的正确分离。在大多数真核生物（如酵母、植物、哺乳动物）中，I 类交叉占总交叉事件的 70%–90%，是减数分裂重组的主导途径。

II 类交叉（Class II crossovers）由不依赖 ZMM 蛋白的替代途径产生，例如 Mus81-Eme1/Mms4 等核酸内切酶介导的机制。这类交叉的发生不存在干扰，即各交叉事件相互独立，在染色体上呈随机分布。尽管所占比例较小，但 II 类交叉作为备用途径，在 I 类交叉功能缺失时仍能维持部分重组，且广泛存在于从酵母到植物的多种生物中，体现了其作为保守重组修复机制的作用。

## Filter haploid cells in scRNA-seq

用于从混合花粉单细胞 RNA-seq 数据中过滤出高质量单倍体细胞，并检测减数分裂交叉互换事件。

文章：[Meiotic recombination dynamics in plants with repeat-based holocentromeres shed light on the primary drivers of crossover patterning](https://www.nature.com/articles/s41477-024-01625-y#Sec1)
GitHub：[detectCO_by_scRNAseq](https://github.com/Raina-M/detectCO_by_scRNAseq)

### 构建参考 SNP 标记集

使用 minimap2 将 Hi-Fi reads 比对到单倍型 hap2，使用 bcftools 进行 call SNP。

```bash
samtools faidx /home/jiaoyuan/04.ref/Palv.hap2.chr.fa

(samtools fastq /home/jiaoyuan/01.projects/snp/Palv_hifi-1.ccs.bam; samtools fastq /home/jiaoyuan/01.projects/snp/Palv_hifi-2.ccs.bam) | \
minimap2 -t 8 -ax map-hifi -R '@RG\tID:Palv_hap2\tSM:Palv\tPL:PACBIO' /home/jiaoyuan/04.ref/Palv.hap2.chr.fa - | \
samtools sort -@ 4 -m 4G -o /home/jiaoyuan/01.projects/snp/Palv_hap2_aligned_sorted.bam

samtools index -@ 4 /home/jiaoyuan/01.projects/snp/Palv_hap2_aligned_sorted.bam

bcftools mpileup -B --threads 32 -Ou -q 20 -a AD,DP -f /public/home/jiaoyuan/01.ref/Palv.hap2.chr.fa /public/home/jiaoyuan/02.crossover/Palv_hap2_aligned_sorted.bam | bcftools call --threads 32 -mv -Oz -o /public/home/jiaoyuan/02.crossover/Palv.snps.hap2.vcf.gz

bcftools index --threads 64 /public/home/jiaoyuan/02.crossover/Palv.snps.hap2.vcf.gz
```

### 单细胞数据处理

**BAM 文件拆分**

花粉单细胞数据使用 DNBC4tools 比对与定量，输出包含细胞标签（CB:Z:）的 BAM 文件。从该 BAM 文件中统计每个细胞 barcode 的 reads 数，保留 reads > 5,000 的有效细胞，并将每个细胞拆分为独立的 BAM 文件，用于后续单细胞分析。

```bash
samtools view anno_decon_sorted.bam | grep -oP 'CB:Z:\S+' | sort | uniq -c | awk '{print $2, $1}' > cell_reads.txt
awk '$2>5000 {print $1}' cell_reads.txt | sed 's/CB:Z://' > good_cells.list
samtools view -H anno_decon_sorted.bam > header.sam
mkdir -p demultiplexed_cells
while read cid; do
    samtools view anno_decon_sorted.bam | awk -v cid="CB:Z:$cid" -v header="header.sam" 'BEGIN{while(getline < header) print} $0~cid{print}' | samtools view -bS - > demultiplexed_cells/${cid}.bam
done < good_cells.list
rm header.sam
```

**去除 PCR 重复**

使用 UMIcollapse 基于 UMI 去除每个细胞 BAM 文件中的 PCR 重复。

```bash
mkdir -p /public/home/jiaoyuan/02.crossover/Pal_pollen_1/dedup_bam
cd /public/home/jiaoyuan/02.crossover/Pal_pollen_1/demux_bam
for bam in *.bam; do
    cid=${bam%.bam}
    umicollapse bam -i ${bam} -o ../dedup_bam/${cid}.dedup.bam
    samtools index ../dedup_bam/${cid}.dedup.bam
done
```

**单细胞 SNP calling**

成熟花粉粒是雄配子体，是单倍体，基因组只包含一套染色体。单细胞 RNA-seq 获得的转录本理论上应来自单一亲本单倍型，除非发生交叉互换。

使用 bcftools call 以单倍体模式（ploidy=1）检测变异，提取 DP4 字段计算参考等位基因和替代等位基因的 reads 数，使用 `get_subset.pl` 将细胞 SNP 与参考 SNP 标记集取交集，获得每个细胞在标记位点上的等位基因计数，所得数据将用于后续单倍型溯源分析。

```bash
bcftools mpileup -Oz -o ${BC}.mpileup.gz -f reference_hap1.fa ${BC}.sorted.dedup.bam

echo -e "*\t*\t*\t*\t1" > ploidy.txt

bcftools call -Am -Ov --ploidy-file ploidy.txt ${BC}.mpileup.gz > ${BC}.vcf

bcftools query -f '%CHROM %POS %REF %ALT %QUAL [ %INDEL %DP %DP4]\n' ${BC}.vcf -o ${BC}.comma.txt

tr ',' '\t' < ${BC}.comma.txt > ${BC}.tabbed.txt

awk '{if ($1 !="ChrC" && $1!="ChrM") print $1 "\t" $2 "\t" $3 "\t" $8+$9 "\t"  $4 "\t"  $10+$11}' \
    ${BC}.tabbed.txt > ${BC}.input

perl get_subset.pl ${BC}.input 1,2 reference_markers.txt 2,3 0 > ${BC}_input_corrected.txt
```

**筛选单倍型细胞**

对每个细胞计算有效标记数（ref+alt>0 的位点数），根据等位基因频率（ref/(ref+alt)）初步判定每个标记的基因型（<0.2 为 ref，>0.8 为 alt，其余为不确定），统计相邻标记基因型发生变化的次数并计算切换率（基因型切换次数 / 总标记数），剔除标记数<400 或切换率>0.07 的双细胞（doublets），所得高质量细胞将用于后续单倍型溯源分析。

```bash
while read BC
do
    awk '$4!=0 || $6!=0' ${BC}_input_corrected.txt > ${BC}_input.tmp
    marker_num=`wc -l ${BC}_input.tmp | cut -d" " -f1`
    awk '{
        if      ( $4 / ($4+$6) < 0.2) print 0;
        else if ( $4 / ($4+$6) > 0.8) print 1;
        else print 0.5;
    }' ${BC}_input.tmp > ${BC}_smt_genotypes.txt
    head -n-1 ${BC}_smt_genotypes.txt > tmp1
    tail -n+2 ${BC}_smt_genotypes.txt > tmp2
    switch_num=`paste tmp1 tmp2 | awk '$2-$1!=0' | wc -l`
    echo -e $BC"\t"$marker_num"\t"$switch_num >> switches.stats 
    rm ${BC}_input.tmp
    rm ${BC}_smt_genotypes.txt
    rm tmp*
done < barcode.list

awk '$2>=400' switches.stats | awk '$3/$2<=0.07 {print $1}' > barcode_gt400markers_no_doublets.list
```

## Crossover Calling

单倍体配子中，如果某条染色体的部分区域来自 hap1、部分来自 hap2，表明在这个花粉形成前的减数分裂中，该染色体发生了交叉互换。通过 scRNA-seq 检测到的一条染色体上单倍型来源的两次转换，即可推断一个交叉互换事件。对每个染色体的标记进行等位基因频率平滑（前后各 2 个标记），第二轮平滑将窗口中多数基因型作为校正依据，根据平滑后的基因型（0/1/0.5）将连续相同基因型的标记合并为基因型块，过滤长度≥1 Mb 且标记数≥5 的块，识别相邻不同基因型块之间的边界并在边界区间内精细定位 CO 断点，输出 CO 位置文件及基因型块文件，所得交叉互换事件将用于后续遗传图谱构建。

使用 `hapCO_identification.R` 批量运行：

```bash
while read BC
do
    Rscript hapCO_identification.R -i ${BC}_input_corrected.txt \
                                   -p $BC \
                                   -c 400 \
                                   -g reference_hap1.genome \
                                   -o outdir
done < barcode_gt400markers_no_doublets.list
```

## 新疆杨花粉细胞交叉互换的识别

### 00.HiFi reads 比对到参考基因组

```bash
samtools faidx /home/jiaoyuan/04.ref/Palv.hap2.chr.fa  
  
(samtools fastq /home/jiaoyuan/01.projects/snp/Palv_hifi-1.ccs.bam; samtools fastq /home/jiaoyuan/01.projects/  
snp/Palv_hifi-2.ccs.bam) | \  
minimap2 -t 8 -ax map-hifi -R '@RG\tID:Palv_hap2\tSM:Palv\tPL:PACBIO' /home/jiaoyuan/04.ref/Palv.hap2.chr.fa -  
| \  
samtools sort -@ 4 -m 4G -o /home/jiaoyuan/01.projects/snp/Palv_hap2_aligned_sorted.bam
samtools index -@ 4 /home/jiaoyuan/01.projects/snp/Palv_hap2_aligned_sorted.bam  
```

### 01.按照染色体 call SNP

```bash
head -19 //home/jiaoyuan/04.ref/Palv.hap2.chr.fa.fai | cut -f1 > /home/jiaoyuan/01.projects/covercross/chr_lis  
t.txt  
grep "^UN" //home/jiaoyuan/04.ref/Palv.hap2.chr.fa.fai | cut -f1 > /home/jiaoyuan/01.projects/covercross/un_li  
st.txt  
cat /home/jiaoyuan/01.projects/covercross/chr_list.txt > /home/jiaoyuan/01.projects/covercross/regions.txt  
un_region=$(cat /home/jiaoyuan/01.projects/covercross/un_list.txt | tr '\n' ',' | sed 's/,$//')  
echo $un_region >> /home/jiaoyuan/01.projects/covercross/regions.txt
```

```bash
mkdir -p /home/jiaoyuan/01.projects/covercross/vcf_split  
mkdir -p /home/jiaoyuan/01.projects/covercross/logs  
export TMPDIR=/home/jiaoyuan/tmp  
mkdir -p /home/jiaoyuan/tmp  
  
samtools flagstat --threads 20 /home/jiaoyuan/01.projects/covercross/Palv_hap2_aligned_sorted.bam > /home/jiao  
yuan/01.projects/covercross/logs/mapping_stats.txt  
samtools coverage /home/jiaoyuan/01.projects/covercross/Palv_hap2_aligned_sorted.bam -o /home/jiaoyuan/01.proj  
ects/covercross/logs/coverage_depth.txt  
  
seq 20 | xargs -I {} -P 20 bash -c '  
   REGION=$(sed -n "{}p" /home/jiaoyuan/01.projects/covercross/regions.txt)  
   bcftools mpileup -Ou -q 10 -Q 10 -a AD,DP -m 2 -F 0.02 -B -r "$REGION" -f /home/jiaoyuan/04.ref/Palv.hap2.  
chr.fa /home/jiaoyuan/01.projects/covercross/Palv_hap2_aligned_sorted.bam | \  
   bcftools call --threads 4 -mv -a GQ -P 0.01 -Oz -o /home/jiaoyuan/01.projects/covercross/vcf_split/part_{}  
.vcf.gz  
   bcftools index -t /home/jiaoyuan/01.projects/covercross/vcf_split/part_{}.vcf.gz  
'  
  
ls /home/jiaoyuan/01.projects/covercross/vcf_split/part_*.vcf.gz | sort -V > /home/jiaoyuan/01.projects/coverc  
ross/vcf_split/merge_list.txt  
  
bcftools concat --threads 20 -f /home/jiaoyuan/01.projects/covercross/vcf_split/merge_list.txt -Ou | \  
bcftools filter -e 'QUAL<10 || INFO/DP<2' -Oz -o /home/jiaoyuan/01.projects/covercross/all.hap2.chr.snp.gz
bcftools index --threads 20 /home/jiaoyuan/01.projects/covercross/all.hap2.chr.snp.gz
bcftools view -i 'GT="0/1" && QUAL>=30' /public/home/jiaoyuan/02.crossover/all.hap2.chr.snp.vcf.gz -Oz -o /public/home/jiaoyuan/02.crossover/f1_hetero_only.vcf.gz
bcftools index -t /public/home/jiaoyuan/02.crossover/f1_hetero_only.vcf.gz
```

### 02.单细胞数据按照细胞进行拆分

```bash
mkdir -p /public/home/jiaoyuan/02.crossover/demux_bam/Pal_pollen_1
cd /public/home/jiaoyuan/02.crossover/demux_bam/Pal_pollen_1

samtools view -@ 32 /public/home/jiaoyuan/03.palv_scrna/dnbc4_pollen/Pal_pollen_1/outs/anno_decon_sorted.bam | awk -F'CB:Z:' 'NF>1 {split($2,a,"[\t ]"); count[a[1]]++} END {for(i in count) if(count[i]>5000) print i}' > list.txt

python -c '
import pysam
valid_cbs = set(line.strip() for line in open("list.txt"))
with pysam.AlignmentFile("/public/home/jiaoyuan/03.palv_scrna/dnbc4_pollen/Pal_pollen_1/outs/anno_decon_sorted.bam", "rb") as fin, pysam.AlignmentFile("filtered.bam", "wb", template=fin) as fout:
    for read in fin:
        if read.has_tag("CB") and read.get_tag("CB") in valid_cbs:
            fout.write(read)
'

samtools sort -@ 32 -t CB filtered.bam -o filtered_sorted.bam
samtools split -@ 32 -d CB -f '%!.bam' -M 50000 filtered_sorted.bam
rm list.txt filtered.bam filtered_sorted.bam

mkdir -p /public/home/jiaoyuan/02.crossover/demux_bam/Pal_pollen_2
cd /public/home/jiaoyuan/02.crossover/demux_bam/Pal_pollen_2

samtools view -@ 32 /public/home/jiaoyuan/03.palv_scrna/dnbc4_pollen/Pal_pollen_2/outs/anno_decon_sorted.bam | awk -F'CB:Z:' 'NF>1 {split($2,a,"[\t ]"); count[a[1]]++} END {for(i in count) if(count[i]>5000) print i}' > list.txt

python -c '
import pysam
valid_cbs = set(line.strip() for line in open("list.txt"))
with pysam.AlignmentFile("/public/home/jiaoyuan/03.palv_scrna/dnbc4_pollen/Pal_pollen_2/outs/anno_decon_sorted.bam", "rb") as fin, pysam.AlignmentFile("filtered.bam", "wb", template=fin) as fout:
    for read in fin:
        if read.has_tag("CB") and read.get_tag("CB") in valid_cbs:
            fout.write(read)
'

samtools sort -@ 32 -t CB filtered.bam -o filtered_sorted.bam
samtools split -@ 32 -d CB -f '%!.bam' -M 50000 filtered_sorted.bam
rm list.txt filtered.bam filtered_sorted.bam
```

### 03.去除每个细胞的 PCR 重复

```bash
ulimit -n 65535

mkdir -p /public/home/jiaoyuan/02.crossover/demux_bam/dedup_bam2/p1
cd /public/home/jiaoyuan/02.crossover/demux_bam/Pal_pollen_1
find . -maxdepth 1 -name "*.bam" | xargs -I {} -P 32 sh -c '
    out_name=$(basename "$1")
    samtools sort -n "$1" | samtools fixmate -m - - | samtools sort - | samtools markdup -r - "/public/home/jiaoyuan/02.crossover/demux_bam/dedup_bam2/p1/$out_name"
    samtools index "/public/home/jiaoyuan/02.crossover/demux_bam/dedup_bam2/p1/$out_name"
' -- {}

mkdir -p /public/home/jiaoyuan/02.crossover/demux_bam/dedup_bam2/p2
cd /public/home/jiaoyuan/02.crossover/demux_bam/Pal_pollen_2
find . -maxdepth 1 -name "*.bam" | xargs -I {} -P 32 sh -c '
    out_name=$(basename "$1")
    samtools sort -n "$1" | samtools fixmate -m - - | samtools sort - | samtools markdup -r - "/public/home/jiaoyuan/02.crossover/demux_bam/dedup_bam2/p2/$out_name"
    samtools index "/public/home/jiaoyuan/02.crossover/demux_bam/dedup_bam2/p2/$out_name"
' -- {}
```

### 04.对每个细胞进行 call SNP 并统计碱基

先将 bam 文件分成每 1000 个放一个文件夹，防止超出 Linux 句柄：

```bash
cd /public/home/jiaoyuan/02.crossover/demux_bam/dedup_bam/p1
ls CELL*.bam | split -l 1000 -d - "part_"
for f in part_*; do
    dir_name=${f#part_}
    mkdir -p "$dir_name"
    cat "$f" | xargs -I {} mv {} {}.bai "$dir_name"/
done
rm part_*

cd /public/home/jiaoyuan/02.crossover/demux_bam/dedup_bam/p2
ls CELL*.bam | split -l 1000 -d - "part_"
for f in part_*; do
    dir_name=${f#part_}
    mkdir -p "$dir_name"
    cat "$f" | xargs -I {} mv {} {}.bai "$dir_name"/
done
rm part_*
```

为每个文件夹生成文件列表索引：

```bash
for d in /public/home/jiaoyuan/02.crossover/demux_bam/dedup_bam/p1/[0-9][0-9]; do ls $d/CELL*.bam > $d/bam.list; sed 's/\.bam$//' $d/bam.list | xargs -n 1 basename > $d/sample.id; done

for d in /public/home/jiaoyuan/02.crossover/demux_bam/dedup_bam/p2/[0-9][0-9]; do ls $d/CELL*.bam > $d/bam.list; sed 's/\.bam$//' $d/bam.list | xargs -n 1 basename > $d/sample.id; done
```

对每个细胞 call SNP

```bash
ulimit -n 100000

for dir in /public/home/jiaoyuan/02.crossover/demux_bam/dedup_bam/p1/[0-9][0-9]; do
    mkdir -p /public/home/jiaoyuan/02.crossover/demux_bam/cellsnp_results/p1/$(basename "$dir")
    cellsnp-lite -S "$dir/bam.list" \
                 -i "$dir/sample.id" \
                 -T /public/home/jiaoyuan/02.crossover/all.hap2.chr.snp.vcf.gz \
                 -O /public/home/jiaoyuan/02.crossover/demux_bam/cellsnp_results/p1/$(basename "$dir") \
                 -p 32 --minCOUNT 1 --minMAF 0 --cellTAG None --UMItag None --gzip --genotype
done

for dir in /public/home/jiaoyuan/02.crossover/demux_bam/dedup_bam/p2/[0-9][0-9]; do
    mkdir -p /public/home/jiaoyuan/02.crossover/demux_bam/cellsnp_results/p2/$(basename "$dir")
    cellsnp-lite -S "$dir/bam.list" \
                 -i "$dir/sample.id" \
                 -T /public/home/jiaoyuan/02.crossover/all.hap2.chr.snp.vcf.gz \
                 -O /public/home/jiaoyuan/02.crossover/demux_bam/cellsnp_results/p2/$(basename "$dir") \
                 -p 32 --minCOUNT 1 --minMAF 0 --cellTAG None --UMItag None --gzip --genotype
done
```

### 05.比对单个细胞的 SNP 与参考 SNP 筛选单倍体细胞

```bash
mkdir -p /public/home/jiaoyuan/02.crossover/high_quality_single_cell/p1/bamfiles
mkdir -p /public/home/jiaoyuan/02.crossover/high_quality_single_cell/p2/bamfiles
```

```python
import os
import numpy as np
import pandas as pd
from scipy.io import mmread

def run_filter(base_dir, bam_base, output_stats, hq_list):
    batches = sorted([d for d in os.listdir(base_dir) if d.isdigit()])
    with open(output_stats, "w") as f_out, open(hq_list, "w") as f_hq:
        for b in batches:
            b_path = os.path.join(base_dir, b)
            ad_path = os.path.join(b_path, "cellSNP.tag.AD.mtx")
            dp_path = os.path.join(b_path, "cellSNP.tag.DP.mtx")
            sample_path = os.path.join(b_path, "cellSNP.samples.tsv")
            
            if not os.path.exists(ad_path):
                continue
                
            ad = mmread(ad_path).tocsc()
            dp = mmread(dp_path).tocsc()
            samples = pd.read_csv(sample_path, header=None)[0].tolist()

            for i in range(len(samples)):
                bc = samples[i]
                cell_dp = dp[:, i].toarray().flatten()
                cell_ad = ad[:, i].toarray().flatten()
                
                mask = cell_dp > 0
                valid_dp = cell_dp[mask]
                valid_ad = cell_ad[mask]
                marker_num = len(valid_dp)
                
                if marker_num > 0:
                    ref_counts = valid_dp - valid_ad
                    ratios = ref_counts / valid_dp
                    genotypes = np.where(ratios < 0.2, 0, np.where(ratios > 0.8, 1, 0.5))
                    switches = np.sum(np.diff(genotypes) != 0) if marker_num > 1 else 0
                else:
                    switches = 0
                
                f_out.write(f"{bc}\t{marker_num}\t{switches}\n")
                
                if marker_num >= 400 and (switches / marker_num <= 0.07):
                    src_bam = os.path.join(bam_base, b, f"{bc}.bam")
                    f_hq.write(f"{src_bam}\n")

run_filter("/public/home/jiaoyuan/02.crossover/demux_bam/cellsnp_results/p1", 
           "/public/home/jiaoyuan/02.crossover/demux_bam/dedup_bam/p1", 
           "/public/home/jiaoyuan/02.crossover/demux_bam/cellsnp_results/p1/switches.stats", 
           "/public/home/jiaoyuan/02.crossover/demux_bam/cellsnp_results/p1/hqsc_bam.list")

run_filter("/public/home/jiaoyuan/02.crossover/demux_bam/cellsnp_results/p2", 
           "/public/home/jiaoyuan/02.crossover/demux_bam/dedup_bam/p2", 
           "/public/home/jiaoyuan/02.crossover/demux_bam/cellsnp_results/p2/switches.stats", 
           "/public/home/jiaoyuan/02.crossover/demux_bam/cellsnp_results/p2/hqsc_bam.list")
```

```bash
cat /public/home/jiaoyuan/02.crossover/demux_bam/cellsnp_results/p1/hqsc_bam.list | xargs -I {} cp {} {}.bai /public/home/jiaoyuan/02.crossover/high_quality_single_cell/p1/bamfiles/
cat /public/home/jiaoyuan/02.crossover/demux_bam/cellsnp_results/p2/hqsc_bam.list | xargs -I {} cp {} {}.bai /public/home/jiaoyuan/02.crossover/high_quality_single_cell/p2/bamfiles/
```

### 06.识别交叉互换

```bash
mkdir -p /public/home/jiaoyuan/02.crossover/high_quality_single_cell/p1/co_results
mkdir -p /public/home/jiaoyuan/02.crossover/high_quality_single_cell/p2/co_results
```

```python
import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
from scipy.io import mmread

def smoother(df):
    df = df.copy()
    df['rawAF'] = df['V4'] / (df['V4'] + df['V6'])
    df['smtAF'] = df['rawAF'].copy()
    window_af = df['rawAF'].rolling(window=5, center=True, min_periods=1).mean()
    df.loc[window_af > 0.7, 'smtAF'] = 1
    df.loc[window_af < 0.3, 'smtAF'] = 0
    df['genotype'] = np.where(df['smtAF'] > 0.8, 1, np.where(df['smtAF'] < 0.2, 0, 0.5))
    return df

def get_blocks(df):
    if df.empty: return pd.DataFrame()
    df = df.copy()
    df['change'] = df['genotype'].diff().ne(0)
    df['block_id'] = df['change'].cumsum()
    blks = []
    for _, g in df.groupby('block_id'):
        blks.append([g.iloc[0]['V2'], g.iloc[-1]['V2'], g.iloc[0]['genotype'], len(g)])
    res = pd.DataFrame(blks, columns=['sta', 'end', 'genotype', 'markernum'])
    return res[(res['end'] - res['sta'] >= 1000000) & (res['markernum'] >= 5)]

def get_bps(blocks, df, chrom):
    bps = []
    if len(blocks) < 2: return pd.DataFrame()
    for i in range(len(blocks) - 1):
        r1, r2 = blocks.iloc[i], blocks.iloc[i+1]
        if r1['genotype'] != r2['genotype'] and r1['genotype'] != 0.5 and r2['genotype'] != 0.5:
            inter = df[(df['V2'] >= r1['sta']) & (df['V2'] <= r2['end'])]
            if len(inter) < 2: continue
            seq = inter['smtAF'].values
            best_s, best_p = -1, 0
            for p in range(len(seq)-1):
                l0, l1 = np.sum(seq[:p+1]==0), np.sum(seq[:p+1]==1)
                r0, r1 = np.sum(seq[p+1:]==0), np.sum(seq[p+1:]==1)
                s = max((l0/(p+1))*(r1/(len(seq)-p-1)), (l1/(p+1))*(r0/(len(seq)-p-1)))
                if s > best_s: best_s, best_p = s, p
            bps.append([chrom, inter.iloc[best_p]['V2'], inter.iloc[best_p+1]['V2'], r1['genotype'], r2['genotype']])
    return pd.DataFrame(bps, columns=['chr', 'sta', 'end', 'pre_gt', 'nxt_gt'])

def plot_co(df, blocks, chr_dict, bc, out):
    chrs = sorted(chr_dict.keys())
    fig, axes = plt.subplots(len(chrs), 1, figsize=(12, 2*len(chrs)), sharex=True)
    if len(chrs) == 1: axes = [axes]
    for i, c in enumerate(chrs):
        ax = axes[i]
        sub = df[df['V1'] == c]
        ax.vlines(sub[sub['V4']>0]['V2'], 0, sub[sub['V4']>0]['V4'], colors='red', linewidth=0.5)
        ax.vlines(sub[sub['V6']>0]['V2'], 0, -sub[sub['V6']>0]['V6'], colors='blue', linewidth=0.5)
        if not blocks.empty:
            c_blks = blocks[blocks['chr']==c]
            for _, b in c_blks.iterrows():
                ax.add_patch(Rectangle((b['sta'], 26), b['end']-b['sta'], 4, color='red' if b['genotype']==1 else 'blue' if b['genotype']==0 else 'purple'))
        ax.set_ylim(-30, 35)
        ax.set_ylabel(c, rotation=0, labelpad=20)
    plt.tight_layout()
    plt.savefig(out)
    plt.close()

def identify_crossover(cellsnp_res_dir, hq_list_file, fai_file, out_dir):
    chr_dict = pd.read_csv(fai_file, sep="\t", header=None, index_col=0)[1].to_dict()
    hq_set = set(pd.read_csv(hq_list_file, header=None)[0].apply(lambda x: os.path.basename(x).replace(".bam","")))
    
    for b in sorted([d for d in os.listdir(cellsnp_res_dir) if d.isdigit()]):
        b_path = os.path.join(cellsnp_res_dir, b)
        vcf_info = pd.read_csv(os.path.join(b_path, "cellSNP.base.vcf.gz"), sep="\t", comment="#", header=None, usecols=[0,1,3,4])
        vcf_info.columns = ['V1', 'V2', 'V3', 'V5']
        samps = pd.read_csv(os.path.join(b_path, "cellSNP.samples.tsv"), header=None)[0].tolist()
        ad_mtx = mmread(os.path.join(b_path, "cellSNP.tag.AD.mtx")).tocsc()
        dp_mtx = mmread(os.path.join(b_path, "cellSNP.tag.DP.mtx")).tocsc()
        
        for i, bc in enumerate(samps):
            if bc not in hq_set: continue
            cad, cdp = ad_mtx[:, i].toarray().flatten(), dp_mtx[:, i].toarray().flatten()
            mask = cdp > 0
            cell_df = vcf_info[mask].copy()
            cell_df['V4'] = cdp[mask] - cad[mask]
            cell_df['V6'] = cad[mask]
            
            blist, bplist = [], []
            for chrom in cell_df['V1'].unique():
                sdf = smoother(cell_df[cell_df['V1']==chrom])
                blks = get_blocks(sdf)
                if not blks.empty:
                    bps = get_bps(blks, sdf, chrom)
                    blks['chr'] = chrom
                    blist.append(blks)
                    bplist.append(bps)
            
            f_blks = pd.concat(blist) if blist else pd.DataFrame()
            f_bps = pd.concat(bplist) if bplist else pd.DataFrame()
            f_bps.to_csv(f"{out_dir}/{bc}_co_pred.txt", sep="\t", index=False)
            plot_co(cell_df, f_blks, chr_dict, bc, f"{out_dir}/{bc}_co_plot.pdf")

identify_crossover("/public/home/jiaoyuan/02.crossover/demux_bam/cellsnp_results/p1", 
                   "/public/home/jiaoyuan/02.crossover/demux_bam/cellsnp_results/p1/hqsc_bam.list", 
                   "/public/home/jiaoyuan/01.ref/Palv.hap2.chr.fa.fai", 
                   "/public/home/jiaoyuan/02.crossover/high_quality_single_cell/p1/co_results")

identify_crossover("/public/home/jiaoyuan/02.crossover/demux_bam/cellsnp_results/p2", 
                   "/public/home/jiaoyuan/02.crossover/demux_bam/cellsnp_results/p2/hqsc_bam.list", 
                   "/public/home/jiaoyuan/01.ref/Palv.hap2.chr.fa.fai", 
                   "/public/home/jiaoyuan/02.crossover/high_quality_single_cell/p2/co_results")
```