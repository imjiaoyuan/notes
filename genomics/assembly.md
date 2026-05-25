# Genome Assembly

## Background

基因组组装（Genome Assembly）是生物信息学领域的核心问题，其本质是将测序仪产生的读取片段（Reads）经过拼接组装，还原成物种完整的染色体序列。基因组组装对于后续分析至关重要，随着测序技术的迭代，组装策略也在发生巨变：早期的二代测序（NGS）序列较短，难以跨越重复区域；而目前主流的三代测序技术（TGS，如 PacBio HiFi 和 Oxford Nanopore）能够获得长达 15kb 甚至 100kb 以上的超长读长，极大地降低了组装难度，提高了基因组的连续性和完整性。
宏观来说，基因组组装可以分为从头组装（De novo assembly）和基于参考基因组的组装（Reference-guided assembly）。从头组装是指不依赖任何已知参考序列，仅凭测序数据构建基因组，这对于解析新物种或检测复杂的结构变异（SV）是必须的；基于参考的组装则是将测序序列比对到近缘物种的基因组骨架上，适合快速重测序分析，但容易丢失物种特有的序列信息。
由于基因组中大量重复序列的存在和生物体倍性的复杂性，组装结果往往难以直接得到完整的染色体。因此，评价组装质量的标准主要涵盖连续性、完整性和准确性：连续性通过 Contig N50 等指标衡量，追求序列片段尽可能长；完整性通常使用 BUSCO 评估基因覆盖度；准确性则要求极低的碱基错误率（高 QV 值）。此外，当前的组装技术已向单倍型分型（Haplotype-resolved）迈进，这是一个更高阶的要求，旨在利用 Hi-C 或长读长数据精准区分父本和母本的两套染色体，甚至追求实现端粒到端粒（T2T）的无缝组装。
在算法层面，基因组组装主要分为基于 De Brujin Graph (DBG) 和基于 OLC (Overlap-Layout-Consensus) 两条路径。DBG 算法将序列打断为 k-mer 构建图谱，运算速度快，适合早期的二代短读长数据，但难以解决长重复序列导致的碎片化问题。相比之下，OLC 及其变体 String Graph 算法直接寻找长序列间的重叠关系，虽然过去因计算量大受限，但随着 PacBio HiFi 等高准确度长读长数据的出现，现代主流软件（如 Hifiasm, Canu）均回归并改良了 OLC 算法，利用长读长跨越重复区并结合 Hi-C 纠错，高效实现高质量的组装。

## Definition

- k-mer：指将一条序列拆解为长度为 k 个碱基的所有连续子字符串。若 read 长度为 L，k-mer 长度为 k，则可产生的 k-mers 数目为 L-k+1。例如序列 AACTGACT，设置 k 为 3，可分割为 AAC、ACT、CTG、TGA、GAC、ACT 共 6 个 k-mers。在基因组拼接中，k 通常设为奇数，以避免在构建 De Bruijn 图时，双链对称的回文序列产生与其自身互补配对的 k-mer，从而导致路径歧义。
- Contig：即片段重叠群。指拼接软件利用 reads 之间的重叠区（overlap）进行组装，获得的连续不含 N 的较长序列。Contig 的长度和准确性是评价基因组初级拼接质量的重要指标。
- Scaffold：即片段框架。由多个已知顺序和方向的 Contigs 组成的序列，Contigs 之间通常存在未知的碱基，以 N 表示（称为 Gap）。在 de novo 组装中，通过双端测序（Paired-end）或大片段文库（Mate-pair）提供的跨度信息，可以确定 Contigs 之间的相对位置关系，从而将其连接成 Scaffold。
- N50：评价基因组拼接质量的加权平均长度指标。将所有 Contigs 按长度从大到小排列并依次累加，当累加长度首次达到总长度的 50% 时，当前加入的那条 Contig 的长度即为 Contig N50。N50 越高，说明拼接的连续性越好。该概念常被误解为序列数量排名第 50 名的长度。同理，若累加至总长度的 90% 则称为 N90；此评估标准同样适用于 Scaffolds。
- 测序深度：是指测序得到的总碱基数与待测基因组大小的比值，可以称为 coverage 或 depth。假设一个基因大小为 2M，测序深度为 50×，那么获得的总数据量为 100M。
- k-mer 深度：每一个 k-mer 出现的频数，也可以用 k-mer 丰度指代。
- 覆盖度：指组装获得的序列占目标基因组的比例（即覆盖广度）。由于基因组中存在高 GC 含量、重复序列及复杂结构，组装结果往往无法覆盖所有区域，这些缺失序列的区域称为 Gap。在单基因组组装前，需重点评估：物种基因组大小、杂合度、GC 含量与分布、重复序列比例、遗传图谱可用性及生物学特性（如细菌的革兰氏属性）。参考近缘物种的已发表基因组是获取这些信息的重要途径。基因组大小的评估不仅关系到对组装完整性的判断，还涉及计算资源需求及测序深度的估算。虽然超大基因组（如>10Gb）对计算机内存和计算性能要求极高，组装难度极大，但随着算法优化和硬件提升，目前已可实现。基因组大小可查阅 IMG、NCBI、动物基因组大小数据库等；若数据库未收录，则需通过流式细胞术、k-mer 分析或福尔根染色等实验手段获得。
- 组装算法：指将测序得到的相对较短的片段（reads）重新拼接还原为原始基因组序列的计算模型。组装的核心任务是利用 reads 间的重叠关系或序列特征构建数学图论模型（Graph），并根据算法从中搜索最优路径，以获得连续的 Contig 序列。目前主流的组装算法分为两种：Overlap-Layout-Consensus (OLC) 算法，主要用于测序读长较长、准确率较高的测序数据；De Bruijn Graph (DBG) 算法，主要用于二代高通量测序产生的短片段数据。

## PacBio + Hi-C 联合组装

### HiFi 组装

BAM 转 fq.gz (如果起始数据是 BAM 格式)

```bash
pbindex /public/home/jiaoyuan/04.peu_gnome/cra5199/rawdata/CRR330913.bam
```

Hifiasm 组装

```bash
hifiasm \
    -o hifi_asm \
    -t 24 \
    --h1 rawdata/CRR330931_f1.fq.gz \
    --h2 rawdata/CRR330931_r2.fq.gz \
    rawdata/CRR330925.bam \
    rawdata/CRR330926.bam \
    rawdata/CRR330927.bam \
    rawdata/CRR330928.bam \
    rawdata/CRR330929.bam \
    rawdata/CRR330930.bam

awk '/^S/{print ">"$2"\n"$3}' hifi_asm.hic.hap1.p_ctg.gfa > hifi_asm.hap1.fasta
awk '/^S/{print ">"$2"\n"$3}' hifi_asm.hic.hap2.p_ctg.gfa > hifi_asm.hap2.fasta
```

- `-o`: 输出文件名前缀
- `-t`: 线程数
- `--h1` 和 `--h2` 输入 Hi-C 数据

flye 组装：

```bash
mkdir -p flye_out

flye \
    --pacbio-raw fastq/CRR330925.fastq.gz fastq/CRR330926.fastq.gz fastq/CRR330927.fastq.gz fastq/CRR330928.fastq.gz \
    --out-dir flye_out \
    --threads 24 \
    --genome-size 500m \
    --asm-coverage 40
```

GFA 转 FASTA

```bash
awk '/^S/{print ">"$2;print $3}' test.bp.p_ctg.gfa > test.p_ctg.fa
```

统计 contig 数量

```bash
grep -c "^S" hifiasm/asm.bp.p_ctg.gfa  
```

提取主要 contigs  

```bash
awk '/^S/{print ">"$2"\n"$3}' hifiasm/asm.bp.p_ctg.gfa > hifiasm/asm.p_ctg.fasta 
``` 
  
提取单倍型 1 和 2 的 contigs  

```bash
awk '/^S/{print ">"$2"\n"$3}' hifi_asm.hic.hap1.p_ctg.gfa > hifi_asm.hap1.fasta
awk '/^S/{print ">"$2"\n"$3}' hifi_asm.hic.hap2.p_ctg.gfa > hifi_asm.hap2.fasta
```

#### 组装质量评估

QUAST 评估

```bash
python /software/quast-5.2.0/quast.py test.fa
```

BUSCO 评估

```bash
## 配置环境变量
export BUSCO_CONFIG_FILE="./config.ini" && \
export AUGUSTUS_CONFIG_PATH="/w/00/g/g05/user681/apps/augustus_config/" && \
export PATH=/data/apps/augustus/augustus-3.2.3/bin:/dta/apps/augustus/augustus-3.2.3:$PATH && \
export PYTHONPATH=/data/apps/BUSCO/BUSCO_3.0.2b/lib/python:$PYTHONPATH

## 运行 BUSCO
python /data/apps/BUSCO/BUSCO_3.0.2b/scripts/run_BUSCO.py \
  -i ../../annotation/dym.hap1/4.gene_convert_v2/Dtu.hap1.chr.v2.longest.pep.fa \
  -m proteins \
  -l /w/00/g/g05/WangDeyan/software/embryophyta_odb10/ \
  -o Dtu.hap1 -c 10
```

- `-i`: 输入的 fasta 文件  
- `-l`: 对应的物种库
- `-o`: 输出目录名
- `-m`: 模式 (proteins/genome/transcriptome)

#### 比对率及覆盖度计算

Minimap2 比对与排序

```bash
minimap2 -ax map-hifi -t 20 ../hap1_ref/Dtu.hap1.chr.v2.fa Dtu.cell1.hifi_reads.fq.gz | \
samtools view -@ 16 -bS | \
samtools sort -@ 16 -o Dtu.hap1.hifi.sorted.bam
```

- `.bam`: 二进制序列比对格式，存储 reads 比对到参考基因组的结果

计算比对率

```bash
samtools flagstat {bam} > {bam}.map
```

计算覆盖率

```bash
samtools depth -aa genome.sort.bam > genome.depth

python /w/00/g/g05/WangYubo/software/genome_check_scripts/stat_coverage.py \
  -i genome.depth \
  -d 1,5,10,20 \
  -o genome.coverage.tsv
```

### Hi-C 辅助组装

#### Juicer 流程

准备运行环境

```bash
conda create -n juicer -c bioconda -c conda-forge bwa samtools python3
```

准备目录结构

```bash
mkdir -p /work/home/aczuuitcbj/jiaoyuan/CRA5200/juicer_work/references
mkdir -p /work/home/aczuuitcbj/jiaoyuan/CRA5200/juicer_work/restrictionsites
mkdir -p /work/home/aczuuitcbj/jiaoyuan/CRA5200/juicer_work/fastq

cp /work/home/aczuuitcbj/jiaoyuan/CRA5200/assembly.fasta /work/home/aczuuitcbj/jiaoyuan/CRA5200/juicer_work/references/
cp /work/home/aczuuitcbj/jiaoyuan/CRA5200/assembly.fasta /work/home/aczuuitcbj/jiaoyuan/CRA5200/juicer_work/restrictionsites/

ln -s /work/home/aczuuitcbj/jiaoyuan/CRA5200/CRR330931_f1.fq.gz /work/home/aczuuitcbj/jiaoyuan/CRA5200/juicer_work/fastq/A_S0_L001_R1_001.fastq.gz
ln -s /work/home/aczuuitcbj/jiaoyuan/CRA5200/CRR330931_r2.fq.gz /work/home/aczuuitcbj/jiaoyuan/CRA5200/juicer_work/fastq/A_S0_L001_R2_001.fastq.gz
ln -sf /work/home/aczuuitcbj/jiaoyuan/juicer-1.6/CPU /work/home/aczuuitcbj/jiaoyuan/juicer-1.6/scripts
```

构建 BWA 索引

```bash
bwa index /work/home/aczuuitcbj/jiaoyuan/CRA5200/juicer_work/references/assembly.fasta
```

生成酶切位点文件

```bash
cd /public/share/ac4w0a7em6/jiaoyuan/CRA5200/juicer_work/restrictionsites/
python3 /work/home/aczuuitcbj/jiaoyuan/juicer-1.6/misc/generate_site_positions.py MboI assembly /work/home/aczuuitcbj/jiaoyuan/CRA5200/juicer_work/restrictionsites/assembly.fasta
```

- 这里 gnome 的名称要和 fasta 文件的名称一致，并且文件会直接生成在 fasta 文件所在的位置

生成染色体长度文件

```bash
awk 'BEGIN{OFS="\t"}{print $1, $NF}' /work/home/aczuuitcbj/jiaoyuan/CRA5200/juicer_work/restrictionsites/assembly_MboI.txt > /work/home/aczuuitcbj/jiaoyuan/CRA5200/juicer_work/references/assembly.chrom.sizes
```

- `genome.fa` 必须在 `references` 目录下
- `restriction_sites`: 酶切位点文件，记录基因组上特定限制性内切酶（如 MboI）的切割位置
- `chrom.sizes`: 染色体长度文件，记录每条染色体 contig 的长度

执行 Juicer

```bash
cd /work/home/aczuuitcbj/jiaoyuan/CRA5200/juicer_work
/work/home/aczuuitcbj/jiaoyuan/juicer-1.6/CPU/juicer.sh \
    -g assembly \
    -s MboI \
    -y /work/home/aczuuitcbj/jiaoyuan/CRA5200/juicer_work/restrictionsites/assembly_MboI.txt \
    -z /work/home/aczuuitcbj/jiaoyuan/CRA5200/juicer_work/references/assembly.fasta \
    -p /work/home/aczuuitcbj/jiaoyuan/CRA5200/juicer_work/references/assembly.chrom.sizes \
    -d /work/home/aczuuitcbj/jiaoyuan/CRA5200/juicer_work \
    -D /work/home/aczuuitcbj/jiaoyuan/juicer-1.6 \
    -t 32
```

- `-g`: 样本/基因组命名
- `-d`: 工作目录 (TopDir)
- `-s`: 酶切类型 (如 MboI)
- `-y`: 酶切位点文件路径 (必须绝对路径)
- `-z`: 参考基因组序列路径 (必须绝对路径)
- `-p`: 染色体长度文件路径 (必须绝对路径)
- `-D`: juicer 软件安装目录
- 所有路径一定要给绝对路径，否则程序会报错

#### 3d-dna 组装 (Scaffolding)

运行流程

```bash
mkdir -p /public/share/ac4w0a7em6/jiaoyuan/CRA5200/3dDNA
cd /public/share/ac4w0a7em6/jiaoyuan/CRA5200/3dDNA

ln -sf /public/share/ac4w0a7em6/jiaoyuan/CRA5200/assembly.fasta /public/share/ac4w0a7em6/jiaoyuan/cra5200/3dDNA/assembly.fasta
ln -sf /public/share/ac4w0a7em6/jiaoyuan/CRA5200/juicer_work/aligned/merged_nodups.txt /public/share/ac4w0a7em6/jiaoyuan/cra5200/3dDNA/merged_nodups.txt

export JUICERTOOLS=/public/share/ac4w0a7em6/jiaoyuan/3ddna/visualize/juicebox_tools.jar

bash /public/share/ac4w0a7em6/jiaoyuan/3ddna/run-asm-pipeline.sh -r 2 -i 10000 /public/share/ac4w0a7em6/jiaoyuan/CRA5200/3dDNA/assembly.fasta /public/share/ac4w0a7em6/jiaoyuan/cra5200/3dDNA/merged_nodups.txt > /public/share/ac4w0a7em6/jiaoyuan/CRA5200/3dDNA/3dDNA.run.log 2>&1
```

- `-r 2`: 纠错轮数（rounds），指在组装过程中尝试纠正错误连接的次数
- `-i 10000`: 输入图谱的阈值（input map threshold），用于过滤过短的 contigs 或低质量区域
- `merged_nodups.txt`: Juicer 运行的核心产出文件，包含去重后的 Hi-C 接触信息，是 3d-dna 组装的必要输入
- `assembly.fa`: 初始组装得到的 Contig 序列文件（如 Hifiasm 的输出）
