# scRNA-seq

## Background

单细胞测序技术的出现是分子生物学领域的一个关键范式转变，提供了个体细胞分辨率的分子信息，极大地深化了我们对生物系统的理解，特别是对细胞异质性和动态变化的认识。传统的批量 RNA 测序（Bulk RNA-seq）通过平均数百万个细胞的转录信号来生成全局表达谱，但这不可避免地掩盖了关键的细微生物学信息，例如稀有细胞群、瞬时细胞状态以及对刺激的微妙反应。单细胞 RNA 测序（scRNA-seq）旨在解决这一根本限制，使研究人员能够以前所未有的精确度定义细胞身份，并追踪细胞在发育、疾病发生和治疗反应等过程中的动态状态变化。近年来，单细胞组学已突破单一转录组模式，迈向多模态分析，联合解析基因组、转录组、蛋白质组（如 ADT）以及空间组学数据，以构建细胞状态的全面分子图谱。

单细胞测序的优势和潜能 主要体现在其高分辨率和对复杂样本的适用性：例如，单细胞核 RNA 测序（snRNA-seq）克服了传统 scRNA-seq 中组织解离导致的细胞脆弱性和应激诱导的转录偏差，是分析大型、易碎（如神经元或心肌细胞）或冷冻存档样本的理想方法。在数据分析方面，人工智能，特别是 Transformer 模型，因其在处理大规模异构数据集上的通用性，正成为单细胞组学分析的潜在变革者。这些模型通过在多样化数据上进行自监督预训练，能够生成上下文特异性的基因表达表征，并创建批次鲁棒性的细胞嵌入（cell embeddings），从而有效地整合跨研究、组织甚至跨物种的数据，实现细胞类型注释的知识迁移。此外，Transformer 还展现出强大的预测能力，例如进行数据插补、跨模态预测，以及通过模拟基因干扰来识别高潜力的治疗靶点。

然而，单细胞测序也面临着固有的挑战和局限性： 测序数据普遍存在高稀疏性（零值过多，即 dropout 事件）和细胞间测序深度的异质性，这需要复杂的生物信息学流程进行质量控制、归一化和插补以消除非生物学噪声。另一个主要的技术挑战是批次效应，即由于实验批次、设备或操作员等非生物学因素造成的系统性差异，这些差异如果不进行有效校正，可能会掩盖真正的生物学变异。具体到不同技术，snRNA-seq 由于主要捕获核 RNA（包括未剪接的转录本），可能表现出更高的稀疏性，需要专门的计算调整。对于微生物单细胞 RNA 测序（smRNA-seq），由于微生物缺乏细胞核、mRNA 半衰期极短且 rRNA 含量高，导致基因表达具有强烈的瞬态和高噪声特征，因此，分析框架必须摒弃简单的真核细胞分析流程，转而采用如基于泛基因组注释（MscT）或基于 K-mer 频率（mKmer）的专用计算策略。最后，尽管 Transformer 模型取得了巨大成功，但它们是否是处理本质上非序列组学数据的最合适架构，以及它们能否全面超越现有领域方法，仍然是该领域内悬而未决的问题，有待进一步的严格基准评估。

## Mapping

### Cell Ranger

Cell Ranger 是 10×genomics 官方出的数据处理分析工具，它的适配度好，功能强大，操作简单等优点，成为了单细胞测序数据分析的必备工具。

下载地址：https://www.10xgenomics.com/support/software/cell-ranger/downloads

下载软件压缩包并解压：

```bash
curl -o cellranger-10.0.0.tar.gz "https://cf.10xgenomics.com/releases/cell-exp/cellranger-10.0.0.tar.gz?Expires=1774731243&Key-Pair-Id=APKAI7S6A5RYOXBWRPDA&Signature=K41JpQFtwHaM9-tV8p6IFE7kw8ELzp8XsAxJhd9Ke1fndbrjzRILzzZmKpwbj1d6lIvg-KIwv~zi-h6GZ7KOpjNmILksxjxTRkuwulr57PxBdnKYOZFDpmuHUJB7GUfuG1nnD~B5~h2sTvZYb8CDicpPfxvPjTEIHZFBiiKaWutjUhqXDrMoqRjC5Gpho7OsIfInMo~37K0PeG0nsXVDQVBG4e4XlE3~ztyIz1jbmk1hA4PY2R72XOarwGsBCoxo1bipmeP2K~Cz6t2aGW9TQst6VDljAgiaGHoz~eDtRt6fmEejxgWH0jaZE3YoDY4UNHjwqYUhQnfVXAT0Z6DSCw__"
tar -xzvf cellranger-10.0.0.tar.gz
```

然后通过 `export` 将软件目录临时加入环境变量调用：

```bash
export PATH=/public/share/ac4w0a7em6/jiaoyuan/cellranger-10.0.0:$PATH
```

准备参考基因组：

```bash
mkdir -p /public/share/ac4w0a7em6/jiaoyuan/scrna/formatted_fastq

# 将 GFF3 格式的注释文件转换为 GTF 格式，cellranger 必须使用 gtf 格式
gffread /public/share/ac4w0a7em6/jiaoyuan/scrna/ref/Palv.hap2.chr.longest.gff3 -T -o /public/share/ac4w0a7em6/jiaoyuan/scrna/ref/Palv.hap2.chr.longest.gtf
cd /public/share/ac4w0a7em6/jiaoyuan/scrna/ref/
# 使用 cellranger 构建参考基因组索引
cellranger mkref --genome=Pal_ref --fasta=/public/share/ac4w0a7em6/jiaoyuan/scrna/ref/Palv.hap2.chr.fa --genes=/public/share/ac4w0a7em6/jiaoyuan/scrna/ref/Palv.hap2.chr.longest.gtf --nthreads=64
```

为 `fastq.gz` 文件创建软链：

```bash
ln -sf /public/share/ac4w0a7em6/jiaoyuan/scrna/20260327_GZ2026030506-01/Pal_pollen_1-CDNA/E250155743_L01_35_1.fq.gz /public/share/ac4w0a7em6/jiaoyuan/scrna/formatted_fastq/Pal_pollen_1_S1_L001_R1_001.fastq.gz
ln -sf /public/share/ac4w0a7em6/jiaoyuan/scrna/20260327_GZ2026030506-01/Pal_pollen_1-CDNA/E250155743_L01_35_2.fq.gz /public/share/ac4w0a7em6/jiaoyuan/scrna/formatted_fastq/Pal_pollen_1_S1_L001_R2_001.fastq.gz
ln -sf /public/share/ac4w0a7em6/jiaoyuan/scrna/20260327_GZ2026030506-01/Pal_pollen_2-CDNA/E250155743_L01_36_1.fq.gz /public/share/ac4w0a7em6/jiaoyuan/scrna/formatted_fastq/Pal_pollen_2_S1_L001_R1_001.fastq.gz
ln -sf /public/share/ac4w0a7em6/jiaoyuan/scrna/20260327_GZ2026030506-01/Pal_pollen_2-CDNA/E250155743_L01_36_2.fq.gz /public/share/ac4w0a7em6/jiaoyuan/scrna/formatted_fastq/Pal_pollen_2_S1_L001_R2_001.fastq.gz
```

对样本执行定量分析：

```bash
cd /public/share/ac4w0a7em6/jiaoyuan/scrna/

cellranger count --id=Pal_pollen_1 --transcriptome=/public/share/ac4w0a7em6/jiaoyuan/scrna/ref/Pal_ref --fastqs=/public/share/ac4w0a7em6/jiaoyuan/scrna/formatted_fastq --sample=Pal_pollen_1 --expect-cells=10000 --localcores=64 --localmem=120 --include-introns=true
cellranger count --id=Pal_pollen_2 --transcriptome=/public/share/ac4w0a7em6/jiaoyuan/scrna/ref/Pal_ref --fastqs=/public/share/ac4w0a7em6/jiaoyuan/scrna/formatted_fastq --sample=Pal_pollen_2 --expect-cells=10000 --localcores=64 --localmem=120 --include-introns=true
```

- `--id` 唯一输出 ID，该参数决定了结果输出文件夹的名称
- `--transcriptome` 转录组参考基因组路径，即预先构建好的索引文件夹目录
- `--fastqs` FASTQ 原始序列文件所在的文件夹路径
- `--sample` 样本前缀名，用于在 FASTQ 目录中识别并匹配该样本对应的多个数据文件
- `--expect-cells` 预期捕获的细胞数量，该数值用于辅助算法自动识别和过滤细胞
- `--localcores` 限制程序运行期间占用的最大 CPU 核心数（此处为 64 核）
- `--localmem` 限制程序运行期间占用的最大内存（此处为 120 GB）
- `--include-introns` 包含内含子比对结果，开启后会将比对到内含子区域的 Read 也计入定量，常用于单细胞核测序（snRNA-seq）或提高基因检出率

### SeekSoulTools

SeekSoulTools 是寻因生物开发的一套处理单细胞转录组数据的软件，目前该软件包含两个模块:

- `rna` 模块：用于识别细胞标签 barcode，比对定量，得到可用于下游分析的细胞表达矩阵，之后进行细胞聚类和差异分析，该模块不仅支持 SeekOne®系列试剂盒产出数据，还 可通过对 barcode 的描述，支持多种自定义设计结构；
- `fast` 模块：该模块专门针对 SeekOne® DD 单细胞全序列转录组试剂盒产出的数据，用于对数据进行 barcode 提取，双端 reads 比对，定量，以及全序列数据特有的指标统计。

软件安装：

```bash
mkdir seeksoultools.1.2.0
cd seeksoultools.1.2.0
curl -C - -o seeksoultools.1.2.0.tar.gz "https://seekgene-public.oss-cn-beijing.aliyuncs.com/software/seeksoultools/seeksoultools.1.2.0.tar.gz"
tar zxf seeksoultools.1.2.0.tar.gz

source ./bin/activate
./bin/conda-unpack

export PATH=`pwd`/bin:$PATH
```

然后通过 `export` 临时设置环境变量来调用：

```bash
source /public/share/ac4w0a7em6/jiaoyuan/seeksoultools.1.2.0/bin/activate
export PATH=/public/share/ac4w0a7em6/jiaoyuan/seeksoultools.1.2.0/bin:$PATH
```

创建所需目录：

```bash
mkdir -p /public/share/ac4w0a7em6/jiaoyuan/scrna/ref/Pal_ref_star
mkdir -p /public/share/ac4w0a7em6/jiaoyuan/scrna/seeksoul_result
```

使用 STAR 构建参考基因组索引：

```bash
/public/share/ac4w0a7em6/jiaoyuan/seeksoultools.1.2.0/bin/STAR --runMode genomeGenerate --runThreadN 64 --genomeDir /public/share/ac4w0a7em6/jiaoyuan/scrna/ref/Pal_ref_star --genomeFastaFiles /public/share/ac4w0a7em6/jiaoyuan/scrna/ref/Palv.hap2.chr.fa --sjdbGTFfile /public/share/ac4w0a7em6/jiaoyuan/scrna/ref/Palv.hap2.chr.longest.gtf --genomeSAindexNbases 13
```

- `--runMode genomeGenerate` 运行模式，指定为“构建参考基因组索引”模式
- `--runThreadN 64` 指定程序运行使用的 CPU 线程/核心数（此处为 64 核）
- `--genomeDir` 指定构建好的基因组索引文件的输出存放目录
- `--genomeFastaFiles` 输入的参考基因组 FASTA 格式的 DNA 序列文件路径
- `--sjdbGTFfile` 输入的基因组注释 GTF 格式文件路径，用于提取外显子和剪接点信息
- `--genomeSAindexNbases 13` 后缀数组（Suffix Array）索引的碱基长度，根据基因组总大小调整以优化内存和检索速度（基因组越小该数值需越小）

运行 seeksoultools RNA 分析：

```bash
cd /public/share/ac4w0a7em6/jiaoyuan/scrna/seeksoul_result

seeksoultools rna run \
  --fq1 /public/share/ac4w0a7em6/jiaoyuan/scrna/20260327_GZ2026030506-01/Pal_pollen_1-CDNA/E250155743_L01_35_1.fq.gz \
  --fq2 /public/share/ac4w0a7em6/jiaoyuan/scrna/20260327_GZ2026030506-01/Pal_pollen_1-CDNA/E250155743_L01_35_2.fq.gz \
  --samplename Pal_pollen_1 \
  --genomeDir /public/share/ac4w0a7em6/jiaoyuan/scrna/ref/Pal_ref_star \
  --gtf /public/share/ac4w0a7em6/jiaoyuan/scrna/ref/Palv.hap2.chr.longest.gtf \
  --chemistry DDV2 \
  --core 64 \
  --include-introns \
  --star_path /public/share/ac4w0a7em6/jiaoyuan/seeksoultools.1.2.0/bin/STAR

seeksoultools rna run \
  --fq1 /public/share/ac4w0a7em6/jiaoyuan/scrna/20260327_GZ2026030506-01/Pal_pollen_2-CDNA/E250155743_L01_36_1.fq.gz \
  --fq2 /public/share/ac4w0a7em6/jiaoyuan/scrna/20260327_GZ2026030506-01/Pal_pollen_2-CDNA/E250155743_L01_36_2.fq.gz \
  --samplename Pal_pollen_2 \
  --genomeDir /public/share/ac4w0a7em6/jiaoyuan/scrna/ref/Pal_ref_star \
  --gtf /public/share/ac4w0a7em6/jiaoyuan/scrna/ref/Palv.hap2.chr.longest.gtf \
  --chemistry DDV2 \
  --core 64 \
  --include-introns \
  --star_path /public/share/ac4w0a7em6/jiaoyuan/seeksoultools.1.2.0/bin/STAR
```

- `--fq1` & `--fq2` 分别输入 R1（包含 Cell Barcode 和 UMI）和 R2（包含 cDNA 序列）的 FASTQ 文件路径
- `--samplename` 样本名称，作为输出文件夹的前缀和结果标识
- `--genomeDir` 指定 STAR 基因组索引的存放目录
- `--gtf` 指定基因组注释 GTF 文件的路径
- `--chemistry` 测序试剂盒的化学版本，此处为 `DDV2`
- `--core` 设置程序运行允许使用的 CPU 核心数（此处为 64 核）
- `--include-introns` 开启内含子计数功能，将比对到内含子区域的序列也纳入基因定量统计
- `--star_path` 手动指定 STAR 比对软件可执行文件的绝对路径

### 华大 DNBC4tools 

华大的 C4 平台使用 cDNA 文库和 Oligo 文库，需要通过 DNBC4tools 进行定量。

项目地址：https://github.com/MGI-tech-bioinformatics/DNBelab_C_Series_HT_scRNA-analysis-software

下载 DNBC4tools：

```bash
curl -o dnbc4tools-3.0.tar.gz "ftp://ftp2.cngb.org/pub/CNSA/data7/CNP0008672/Single_Cell/CSE0000574/dnbc4tools-3.0.tar.gz"
tar -xzvf dnbc4tools-3.0.tar.gz
```

构建基因组索引：

```bash
cd /public/share/ac4w0a7em6/jiaoyuan/scrna/ref/
/public/share/ac4w0a7em6/jiaoyuan/dnbc4tools3.0/dnbc4tools rna mkref \
    --ingtf /public/share/ac4w0a7em6/jiaoyuan/scrna/ref/Palv.hap2.chr.longest.gtf \
    --fasta /public/share/ac4w0a7em6/jiaoyuan/scrna/ref/Palv.hap2.chr.fa \
    --threads 64 \
    --species Pal_ref
```

执行 RNA 定量：

```bash
mkdir -p /public/share/ac4w0a7em6/jiaoyuan/scrna/dnbc4_results
cd /public/share/ac4w0a7em6/jiaoyuan/scrna/dnbc4_results

/public/share/ac4w0a7em6/jiaoyuan/dnbc4tools3.0/dnbc4tools rna run \
    --cDNAfastq1 /public/share/ac4w0a7em6/jiaoyuan/scrna/20260327_GZ2026030506-01/Pal_pollen_1-CDNA/E250155743_L01_35_1.fq.gz \
    --cDNAfastq2 /public/share/ac4w0a7em6/jiaoyuan/scrna/20260327_GZ2026030506-01/Pal_pollen_1-CDNA/E250155743_L01_35_2.fq.gz \
    --oligofastq1 /public/share/ac4w0a7em6/jiaoyuan/scrna/20260327_GZ2026030506-01/Pal_pollen_1-oligo/E250155743_L01_83_1.fq.gz \
    --oligofastq2 /public/share/ac4w0a7em6/jiaoyuan/scrna/20260327_GZ2026030506-01/Pal_pollen_1-oligo/E250155743_L01_83_2.fq.gz \
    --genomeDir /public/share/ac4w0a7em6/jiaoyuan/scrna/ref/Pal_ref \
    --name Pal_pollen_1 \
    --threads 64

/public/share/ac4w0a7em6/jiaoyuan/dnbc4tools3.0/dnbc4tools rna run \
    --cDNAfastq1 /public/share/ac4w0a7em6/jiaoyuan/scrna/20260327_GZ2026030506-01/Pal_pollen_2-CDNA/E250155743_L01_36_1.fq.gz \
    --cDNAfastq2 /public/share/ac4w0a7em6/jiaoyuan/scrna/20260327_GZ2026030506-01/Pal_pollen_2-CDNA/E250155743_L01_36_2.fq.gz \
    --oligofastq1 /public/share/ac4w0a7em6/jiaoyuan/scrna/20260327_GZ2026030506-01/Pal_pollen_2-oligo/E250155743_L01_84_1.fq.gz \
    --oligofastq2 /public/share/ac4w0a7em6/jiaoyuan/scrna/20260327_GZ2026030506-01/Pal_pollen_2-oligo/E250155743_L01_84_2.fq.gz \
    --genomeDir /public/share/ac4w0a7em6/jiaoyuan/scrna/ref/Pal_ref \
    --name Pal_pollen_2 \
    --threads 64
```

- `rna run` 执行 DNBC4 单细胞转录组定量分析的主流程
- `--cDNAfastq1` & `--cDNAfastq2` 输入 cDNA 库（通常包含 cDNA 序列）的 R1 和 R2 端 FASTQ 文件
- `--oligofastq1` & `--oligofastq2` 输入 Oligo 库（包含 Cell Barcode 和 UMI 信息）的 R1 和 R2 端 FASTQ 文件
- `--genomeDir` 指定由 `dnbc4tools rna mkref` 预先构建好的基因组参考索引目录
- `--name` 设置样本名称，用于命名输出文件夹和结果文件前缀
- `--threads 64` 指定程序运行允许使用的 CPU 线程数（此处为 64 核）

## Seurat

### 参考文献

-   [Single-nucleus transcriptomics revealed auxin-driven mechanisms of wood plasticity to enhance severe drought tolerance in poplar](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-025-03794-1)
-   [Transcriptional landscape of highly lignified poplar stems at single-cell resolution](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-021-02537-2#Sec2)
-   [Single-cell and spatial multiomics identifies heterogeneous xylem development driven by mechanical stress in Populus](<linkinghub.elsevier.com/retrieve/pii/S1534-5807(25)00436-8>)

### 数据结构

#### 数据存储形式

- 原始表达矩阵
    -   barcodes.tsv.gz: 这个文件包含所有被识别到的细胞的唯一条形码（barcode）。每一行对应一个细胞。
    -   features.tsv.gz (或  genes.tsv.gz): 列出了所有检测到的基因（或其他特征，features）。每一行对应一个基因，通常包含基因 ID 和基因名。
    -   matrix.mtx.gz: 稀疏矩阵，记录了每个基因在每个细胞中的表达计数值（通常是 UMI 计数）。“稀疏”意味着文件中只记录非零的表达值，从而大大节省存储空间
-   R 对象 - RDS 文件：存储单个 R 对象的文件格式。在 Seurat 分析中，创建了一个 Seurat 对象，并对其进行了质控、降维、聚类等一系列分析后，这个包含了所有数据和分析结果的复杂对象可以被保存为一整个.rds 文件。一个.rds 文件通常包含一个完整的 Seurat 对象，其中包括了原始表达矩阵、细胞和基因的元数据（metadata）、归一化后的数据、降维结果（如 PCA, UMAP, tSNE）、细胞聚类信息以及差异表达基因分析结果等 - 一个.RData 文件可以保存当前 R 工作环境中的所有对象，而不仅仅是单个对象。

数据转换

```R
library(Seurat)
## 这个文件夹里应该包含 barcodes.tsv.gz, features.tsv.gz, matrix.mtx.gz
data_directory <- "C:/Users/YourName/Documents/scRNA_data/filtered_feature_bc_matrix/"
expression_data <- Read10X(data.dir = data_directory)
## 创建Seurat对象
pbmc_seurat <- CreateSeuratObject(counts = expression_data,
                                  project = "pbmc3k",
                                  min.cells = 3,
                                  min.features = 200)
## 保存为RDS文件
saveRDS(pbmc_seurat, file = "C:/Users/YourName/Documents/scRNA_data/pbmc3k_raw.rds")
```

#### 查看数据结构

加载 Rdata 数据

```R
file_path <- "/home/jiaoyuan/01.projects/pal-scr/CCA_combined.RData"
load(file_path)
```

查看对象

```R
ls()
```

```R
> ls()
[1] "combined_cca" "file_path"
```

查看 Seurat 对象基本信息

```R
library(Seurat)
print(combined_cca)
```

```R
An object of class Seurat
52237 features across 122313 samples within 2 assays
Active assay: integrated (3000 features, 3000 variable features)
 1 other assay present: RNA
 3 dimensional reductions calculated: pca, umap, tsne
```

- 这个对象总共包含了 122,313 个细胞和 52,237 个基因的数据，这些数据被分成了两种类型储存。
- 当前默认使用的是名为“integrated”的数据集，它只包含 3000 个用于分析的高变异基因，并且已经校正过批次效应。
- 对象里还储存着一个名为“RNA”的数据集，它包含了所有基因未经处理的原始表达数据。
- 数据已经完成了三种降维计算（PCA、UMAP 和 tSNE），可以直接用来画图展示细胞分群。

使用 glimpse() 查看 metadata 数据

```R
library(dplyr)
glimpse(combined_cca@meta.data)
```

```R
Rows: 122,313
Columns: 38
$ orig.ident             <chr> "ssu7d2", "ssu7d2", "ssu7d2", "ssu7d2", "ssu7d2…
$ nCount_RNA             <dbl> 3418, 2811, 4139, 1858, 1170, 1710, 5071, 3719,…
$ nFeature_RNA           <int> 2432, 2141, 2767, 1410, 936, 841, 2714, 2629, 6…
$ pANN_0.25_0.04_434     <dbl> 0.24622030, 0.19438445, 0.15550756, 0.28941685,…
$ doublet_info           <chr> "Singlet", "Singlet", "Singlet", "Singlet", "Si…
$ sample                 <chr> "bud", "bud", "bud", "bud", "bud", "bud", "bud"…
$ pANN_0.25_0.04_361     <dbl> NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA,…
$ split                  <chr> NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA,…
$ pANN_0.25_0.04_329     <dbl> NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA,…
$ pANN_0.25_0.07_432     <dbl> NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA,…
$ pANN_0.25_0.005_410    <dbl> NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA,…
$ pANN_0.25_0.005_208    <dbl> NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA,…
$ pANN_0.25_0.06_408     <dbl> NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA,…
$ pANN_0.25_0.06_359     <dbl> NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA,…
$ pANN_0.25_0.06_370     <dbl> NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA,…
$ pANN_0.25_0.02_153     <dbl> NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA,…
$ pANN_0.25_0.03_163     <dbl> NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA,…
$ pANN_0.25_0.02_193     <dbl> NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA,…
$ pANN_0.25_0.01_62      <dbl> NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA,…
$ pANN_0.25_0.005_175    <dbl> NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA,…
$ pANN_0.25_0.01_65      <dbl> NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA,…
$ pANN_0.25_0.02_96      <dbl> NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA,…
$ pANN_0.25_0.005_113    <dbl> NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA,…
$ pANN_0.25_0.04_13      <dbl> NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA,…
$ pANN_0.25_0.005_96     <dbl> NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA,…
$ pANN_0.25_0.005_172    <dbl> NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA,…
$ pANN_0.25_0.005_195    <dbl> NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA,…
$ pANN_0.25_0.005_212    <dbl> NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA,…
$ pANN_0.25_0.005_192    <dbl> NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA,…
$ pANN_0.25_0.005_200    <dbl> NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA,…
$ pANN_0.25_0.005_174    <dbl> NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA,…
$ pANN_0.25_0.11_194     <dbl> NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA,…
$ pANN_0.25_0.3_89       <dbl> NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA,…
$ S.Score                <dbl> -0.010448613, -0.031776084, 0.109424977, 0.1351…
$ G2M.Score              <dbl> -0.005057777, 0.031466948, -0.032960235, 0.0433…
$ Phase                  <chr> "G1", "G2M", "S", "S", "S", "S", "G1", "S", "G1…
$ integrated_snn_res.0.5 <fct> 2, 2, 0, 4, 7, 0, 0, 1, 6, 3, 11, 7, 7, 5, 13, …
$ seurat_clusters        <fct> 2, 2, 0, 4, 7, 0, 0, 1, 6, 3, 11, 7, 7, 5, 13, …
```

-   `$ orig.ident`: 原始样本 ID，代表这个细胞最初来自哪个文件。
-   `$ nCount_RNA`: 总 RNA 分子数，代表在这个细胞里检测到了多少个 RNA 分子。
-   `$ nFeature_RNA`: 总基因数，代表在这个细胞里检测到了多少种不同的基因。
-   `$ pANN_0.25_0.04_434`: 双细胞预测分数，是算法判断细胞是否为“双细胞”时的一个中间计算分数。
-   `$ doublet_info`: 双细胞判断结果，直接告诉您这个细胞是“Singlet”（单细胞）还是“Doublet”（双细胞）。
-   `$ sample`: 样本分组名。
-   `$ pANN_...` : 双细胞预测的附加分数，通常无需关注。
-   `$ split`: 自定义分组，一个备用列，用于根据需要对细胞进行任何额外的标记或分组。
-   `$ S.Score`: S 期得分，表示细胞与细胞周期 S 期（DNA 复制期）的相似程度。
-   `$ G2M.Score`: G2/M 期得分，一个数值，表示细胞与细胞周期 G2/M 期（细胞分裂期）的相似程度。
-   `$ Phase`: 细胞周期阶段，根据上面两个得分，最终判断细胞处于 G1, S, 还是 G2M 期。
-   `$ integrated_snn_res.0.5`: 聚类结果（分辨率 0.5），表示在分辨率参数为 0.5 时，细胞被分到了哪个簇。
-   `$ seurat_clusters`: 最终聚类结果，这是默认使用的细胞分群结果，显示细胞属于哪个簇（Cluster）。

查看数据层 (assays)

```R
names(combined_cca@assays)
```

```R
[1] "RNA"        "integrated"
```

查看降维结果

```R
names(combined_cca@reductions)
```

```R
[1] "pca"  "umap" "tsne"
```

查看 cluster 名称

```R
levels(combined_cca)
```

```R
 [1] "0"  "1"  "2"  "3"  "4"  "5"  "6"  "7"  "8"  "9"  "10" "11" "12" "13" "14"
[16] "15" "16"
```

查看 cluster 里的细胞数量

```R
table(combined_cca$seurat_clusters)
```

```R
    0     1     2     3     4     5     6     7     8     9    10    11    12
34455 17702 15509 14312 14001  7640  5543  2371  2135  2074  1983   976   909
   13    14    15    16
  904   822   697   280
```

查看不同样本在各个 cluster 中的分布

```R
table(combined_cca$seurat_clusters, combined_cca$sample)
```

```R
       bud cytoledon  leaf  root  stem
  0   1315     14216  8243  8023  2658
  1    167      8660   512  2900  5463
  2   1397      8271  4169   786   886
  3    842      9269  1980  1778   443
  4   1244      5700  2726  2561  1770
  5    790      2605  2163   943  1139
  6     48      1951   113  3194   237
  7   1762       239    47   141   182
  8    108      1213   211   487   116
  9      6       643    70   267  1088
  10    83       758   525   135   482
  11   276       156   110    62   372
  12    53        59   247    21   529
  13   127       184   244    80   269
  14    13       226    19   151   413
  15     9       292    30   235   131
  16     8       239    29     3     1
```

### Seurat

#### v5

##### 设置 Seurat 对象

首先设置 Seurat 对象，这里的数据来自 NCBI GEO 数据集 [GSE190649](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE190649)

```R
library(dplyr)
library(Seurat)
library(patchwork)

## 目录下包含 matrix.mtx, features.tsv (或 genes.tsv)和 barcodes.tsv 文件
bud.data <- Read10X(data.dir = "/home/jiaoyuan/01.projects/pal-scr/bud")
```

初始化 Seurat 对象

```R
bud.seurat <- CreateSeuratObject(counts = bud.data, project = "pal-scr-bud", min.cells = 3, min.features = 200)
bud.seurat
```

```R
An object of class Seurat
31214 features across 10693 samples within 1 assay
Active assay: RNA (31214 features, 0 variable features)
```

-   `project = "pal-scr-bud"` Project 名称，可以通过 `Project(bud.seurat)` 查看或者修改
-   `min.cells = 3` 和 `min.features = 200` 是初始的过滤参数

查看指定的五个基因在前 30 个细胞中的数据

```R
bud.data[c("pal-pou04668", "pal-pou04666", "pal-pou04665", "pal-pou04664", "pal-pou04663"), 1:30]
```

```R
5 x 30 sparse Matrix of class "dgCMatrix"
  [[ suppressing 30 column names ‘AAACCCAAGCCAACCC-1’, ‘AAACCCACAAATTGGA-1’, ‘AAACCCACAACCGACC-1’ ... ]]

pal-pou04668 . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
pal-pou04666 1 . . . 2 . . . . . . . . . . . 1 . 1 . 1 . . . . 1 1 . . .
pal-pou04665 . 2 . . . . . . 3 . . . 2 . . . 2 . . . 1 . . . . 1 1 . . .
pal-pou04664 2 . 1 . 1 . . . . 1 . . . . . . . . . . . . . . . . 1 . . .
pal-pou04663 . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
```

##### 标准预处理流程

Seurat 常用 QC 指标有三个

- 每个细胞的独特基因数量
- 细胞内分子总数
- 线粒体基因组读数百分比
    其中，低质量细胞的独特基因少，细胞双联体或多联体的独特基因多；细胞内分子总数和独特基因数量高度相关；低质量或濒死细胞的线粒体污染更高。另外，可通过 `PercentageFeatureSet()` 函数计算线粒体 QC 指标，计算时会把以 “MT-” 开头的基因当作线粒体基因集。

计算线粒体 QC 指标并添加到元数据

```R
## 处理非模式生物时，线粒体基因的前缀通常不是"MT-"，因此需要实际情况进行修改
## 这里处理的是杨树，所以需要自己提供线粒体和叶绿体的列表进行去除
bud.seurat[["percent.mt"]] <- PercentageFeatureSet(bud.seurat, pattern = "^Ptm-")
```

查看前五个细胞的 QC 指标

```R
head(bud.seurat@meta.data, 5)
```

```R
                    orig.ident nCount_RNA nFeature_RNA percent.mt
AAACCCAAGCCAACCC-1 pal-scr-bud       3418         2432          0
AAACCCACAAATTGGA-1 pal-scr-bud       2811         2141          0
AAACCCACAACCGACC-1 pal-scr-bud       4139         2767          0
AAACCCACACTATGTG-1 pal-scr-bud       1858         1410          0
AAACCCACAGTTAGAA-1 pal-scr-bud       4656         3271          0
```

可视化 QC 指标并用于过滤

```R
## 将QC指标绘制成小提琴图
VlnPlot(bud.seurat, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
```

```R
## FeatureScatter用于可视化特征之间的关系
plot1 <- FeatureScatter(bud.seurat, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(bud.seurat, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2
```

### 流程

#### CCA v4

Seurat 最经典的数据整合流程，其核心思想是寻找数据集之间的相关性最大的线性组合。它将两个数据集投影到一个共享的低维空间，在这个空间里，两组数据的相关性最大化。在此空间内，寻找互为最近邻的细胞对，即寻找锚点 (Anchors)。根据锚点周围的拓扑结构，给锚点打分，去除错误的匹配。

```R
bud <- readRDS("bud/bud.RDS")
cytoledon <- readRDS("cytoledon/cytoledon_hypocotyl.RDS")
leaf <- readRDS("leaf/leaf.RDS")
root <- readRDS("root/root.RDS")
stem <- readRDS("stem/stem.RDS")

bud$sample <- 'bud'
cytoledon$sample <- 'cytoledon'
leaf$sample <- 'leaf'
root$sample <- 'root'
stem$sample <- 'stem'
```

载入五个数目组织的数据，让分别给其打上对应的标签。如果一个样本具有重复，其多个重复之间理应没有批次效应，可以直接进行物理合并，使用 `stem <- merge(stem1, y = c(stem2, stem3))` 进行合并。

```R
combined_list_cca <- list(bud, cytoledon, leaf, root, stem)
combined_list_cca <- lapply(X = combined_list_cca, FUN = function(x) {
    x <- NormalizeData(x)
    x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 3000)
})
```

将五个组织的数据放入一个列表 `combined_list_cca`，逐个样本 (X) 进行处理。先进行标准化 `NormalizeData(x)`：

$$
\frac{1}{2} \ln \left( \frac{count}{TotalCounts} \times 10000 + 1 \right)
$$

将测序深度不同的几个样本强制对其在同一个量级。然后进行方差分析 `FindVariableFeatures(x, selection.method = "vst", nfeatures = 3000)`，使用 vst 方法，计算每个基因在单样本内的均值和方差，寻找出 3000 个高变基因。有一部分管家基因在所有的细胞类型中表达都很高很稳定，所以计算单个样本的中的表达量可以防止样本特有的基因被其他样本中的信号影响。这里选取高变基因的方法有：

-   vst：基于局部多项式回归拟合方差与均值的稳定关系，对单细胞数据最稳定，考虑技术噪音，适合大多数情况，为 Seurat 默认使用的方法。
-   mvp：基于均值 - 方差关系，计算速度快，但是对高表达基因的偏向性比较强。
-   dispersion：基于离散度，传统方法，可能受极端值影响较大。

```R
features_cca <- SelectIntegrationFeatures(object.list = combined_list_cca, nfeatures = 3000)
```

寻找整合锚点，在 cca 列表的数据集中，共同选取 3000 个高变基因作为后续整合的特征集，要求选出的基因在多数样本中都是高变的。

```R
cca.anchors <- FindIntegrationAnchors(object.list = combined_list_cca, anchor.features = features_cca)
```

将多个数据集的表达矩阵投影在一个低维共享空间内，寻找互为最近邻的细胞，其互为一对锚点，然后进行打分，去除错误的匹配。使用 `features_cca` 基因集在多个数据集之间寻找锚点。

```R
combined_cca <- IntegrateData(anchorset = cca.anchors)
DefaultAssay(combined_cca) <- "integrated"
```

`IntegrateData` 方法使用高分锚点对关系计算校正向量，使用校正向量对每个数据集的基因表达矩阵进行批次效应校正，生成一个新的、统一的表达矩阵 (新数值 = 原始数值 - 批次特异性偏差)，存储在 integrated 这个数据层，并且将其作为默认的分析对象。

```R
g2m_genes <- read.table('pal_G2M.txt', header = FALSE, stringsAsFactors = FALSE)$V1
s_genes <- read.table('pal_S.txt', header = FALSE, stringsAsFactors = FALSE)$V1

combined_cca <- CellCycleScoring(
  object = combined_cca,
  s.features = s_genes,
  g2m.features = g2m_genes,
  set.ident = FALSE
)
```

基于提供的 S 期基因列表 `s_genes` 和 G2/M 期基因列表 `g2m_genes`，计算每个细胞的 S 期得分和 G2/M 期得分，得到两个分数：S.Score 和 G2M.Score。根据得分差异，为每个细胞分配一个细胞周期阶段标签（Phase）：

- 如果 S 期得分 > G2/M 期得分，标记为 `S` 期
- 如果 G2/M 期得分 > S 期得分，标记为 `G2M` 期
- 否则标记为 `G1` 期

然后将 Phase 信息添加到 combined_cca 对象的 meta.data 中，此举可以回归掉细胞周期效应，避免细胞周期差异干扰细胞类型鉴定，或者用于分析细胞周期状态与生物学过程的关系。

```R
combined_cca <- ScaleData(combined_cca, vars.to.regress = c("S.Score", "G2M.Score"))
```

对于每一个基因，基于上一步 `CellCycleScoring()` 计算得到的 S 期得分 `S.Score` 和 G2/M 期得分 `G2M.Score` 建立线性模型：

$$
\text{Expression}_g = \beta_{1g} \times \text{S.Score} + \beta_{2g} \times \text{G2M.Score} + \text{Residual}_g
$$

然后提取残差（Residual）并进一步进行 Z-score 标准化（均值为 0，方差为 1）。对数据进行标准化缩放时，消除了细胞周期状态带来的技术变异，使后续分析能更专注于生物学差异。

```R
combined_cca <- RunPCA(combined_cca, npcs = 30, verbose = FALSE)
```

对整合后的数据进行主成分分析（PCA），计算前 30 个主成分，用于降维和后续分析。

```R
combined_cca <- RunUMAP(combined_cca, reduction = "pca", dims = 1:20)
```

使用前 20 个主成分运行 UMAP 降维，生成二维 UMAP 坐标用于可视化。

```R
combined_cca <- RunTSNE(combined_cca, reduction = "pca", dims = 1:20)
```

使用前 20 个主成分运行 t-SNE 降维，生成二维 t-SNE 坐标用于可视化。

```R
combined_cca <- FindNeighbors(combined_cca, reduction = "pca", dims = 1:20)
```

基于前 20 个主成分构建 K 近邻图，计算细胞间的相似性关系，为聚类做准备。

```R
combined_cca <- FindClusters(combined_cca, resolution = 0.5)
```

使用 Louvain 算法在近邻图上进行细胞聚类，分辨率参数 0.5 控制聚类粒度（值越大聚类越多越细）。

```R
DefaultAssay(combined_cca) <- "RNA"
cca_markers <- FindAllMarkers(combined_cca, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.58)
write.csv(cca_markers, file.path(out_dir, "cca_markers.csv"), row.names = FALSE)
```

将默认分析对象切换回原始的 RNA 表达数据（"RNA" assay），在未校正的原始数据上寻找差异基因。为 `FindClusters` 的聚类结果寻找每个 cluster 的特异性标记基因：

-   `only.pos = TRUE`：只保留在 cluster 中上调的基因
-   `min.pct = 0.25`：基因至少在 25% 的细胞中表达
-   `logfc.threshold = 0.58`：对数倍变化阈值约 1.5 倍（$2^{0.58}≈1.5$）

#### SCT v4

SCTransform (Regularized Negative Binomial Regression) 数据整合流程，其核心思想是通过正则化负二项回归模型对 scRNA-seq 数据进行建模。不同于传统的对数标准化，SCTransform 能够更好地去除测序深度（Sequencing Depth）的影响，保留生物学异质性，同时避免了伪计数（pseudocounts）引入的偏差。在整合多组数据时，基于 SCT 修正后的皮尔逊残差（Pearson Residuals）来寻找锚点，通常能获得更清晰的批次校正效果。

```R
combined_list_sct <- list(bud, cytoledon, leaf, root, stem)

for(i in 1:length(combined_list_sct)){
    combined_list_sct[[i]] <- CellCycleScoring(combined_list_sct[[i]], s.features = s_genes, g2m.features = g2m_genes, set.ident = FALSE)
    combined_list_sct[[i]] <- SCTransform(
        combined_list_sct[[i]],
        variable.features.n = 3000,
        vars.to.regress = c("S.Score", "G2M.Score"),
        verbose = FALSE
    )
}
```

将五个组织的数据放入列表 `combined_list_sct`，逐个样本进行处理。这里与 CCA 流程不同，SCTransform 流程替代了传统的 NormalizeData、FindVariableFeatures 和 ScaleData 三个步骤。

首先运行 `CellCycleScoring` 计算 S 期和 G2/M 期得分。随后直接运行 `SCTransform`，该函数对每个基因建立负二项回归模型：

$$
\log(E[\text{count}]) = \beta_0 + \beta_1 \log_{10}(\text{TotalCounts}) + \sum \beta_k \text{Covariates}
$$

其中 Covariates 这里指定为 `vars.to.regress = c("S.Score", "G2M.Score")`。通过将细胞周期得分作为协变量纳入回归模型，可以直接在建模过程中消除细胞周期对基因表达的影响（而不是在后续 ScaleData 步骤中处理）。

模型计算得到的皮尔逊残差（Pearson Residuals）作为后续分析的输入：

$$
\text{Residual} = \frac{\text{observed} - \text{expected}}{\sqrt{\text{expected} + \frac{\text{expected}^2}{\theta}}}
$$

这实现了数据的方差稳定变换（VST），有效地去除了由测序深度带来的技术噪音，同时找出了 3000 个高变特征（variable features）。

```R
features_sct <- SelectIntegrationFeatures(object.list = combined_list_sct, nfeatures = 3000)
combined_list_sct <- PrepSCTIntegration(object.list = combined_list_sct, anchor.features = features_sct)
```

`SelectIntegrationFeatures` 在所有对象中根据 SCT 排序选择前 3000 个最具有变异性的基因用于整合。

紧接着必须运行 `PrepSCTIntegration`。这是一个 SCT 整合特有的步骤。由于 SCTransform 默认只计算高变基因的残差，而整合过程可能需要用到特定选出的 `features_sct` 列表中的基因。该函数确保列表中所有对象对于选定的整合特征集都完整计算了皮尔逊残差，保证后续锚点查找的准确性。

```R
sct.anchors <- FindIntegrationAnchors(object.list = combined_list_sct, normalization.method = "SCT", anchor.features = features_sct)
combined_sct <- IntegrateData(anchorset = sct.anchors, normalization.method = "SCT")
```

寻找锚点并整合数据。这里必须显式指定 `normalization.method = "SCT"`。

-   `FindIntegrationAnchors`：在基于 SCT 残差的空间中寻找互为最近邻的细胞（锚点）。
-   `IntegrateData`：利用锚点信息，对 SCT 残差进行校正，生成一个新的 "integrated" assay。这个 assay 包含了校正后的皮尔逊残差，专门用于后续的降维和聚类，而不是直接用于差异分析。

```R
combined_sct <- RunPCA(combined_sct, npcs = 30, verbose = FALSE)
combined_sct <- RunUMAP(combined_sct, reduction = "pca", dims = 1:20)
combined_sct <- RunTSNE(combined_sct, reduction = "pca", dims = 1:20)
combined_sct <- FindNeighbors(combined_sct, reduction = "pca", dims = 1:20)
combined_sct <- FindClusters(combined_sct, resolution = 0.5)
```

基于整合后的数据（Integrated SCT residuals）进行标准的下游分析流程：

-   `RunPCA`：计算主成分。
-   `RunUMAP` / `RunTSNE`：基于前 20 个主成分进行非线性降维，用于可视化。
-   `FindNeighbors`：构建 KNN 图。
-   `FindClusters`：基于 Louvain 算法进行聚类，分辨率设为 0.5。

```R
pdf(file.path(out_dir, "sct_umap.pdf"), width = 12, height = 10)
print(DimPlot(combined_sct, reduction = "umap", group.by = 'sample', label = TRUE, raster = FALSE))
print(DimPlot(combined_sct, reduction = "umap", group.by = "Phase"))
print(DimPlot(combined_sct, reduction = "umap", label = TRUE, raster = FALSE))
print(DimPlot(combined_sct, reduction = "tsne", label = TRUE, raster = FALSE))
dev.off()

pdf(file.path(out_dir, "sct_umap_sample.pdf"), width = 20, height = 6)
print(DimPlot(combined_sct, reduction = "umap", split.by = "sample", label = TRUE))
dev.off()

save(combined_sct, file = "SCT_combined.RData")
```

可视化与保存结果。

- 生成 `sct_umap.pdf`：展示样本来源分布（查看整合效果）、细胞周期分布（查看回归效果）、以及聚类结果。
- 生成 `sct_umap_sample.pdf`：通过 `split.by = "sample"` 将不同样本拆分展示，便于直观对比各组织中细胞类型的存无或丰度差异。
- 最后保存完整的 Seurat 对象。

```R
DefaultAssay(combined_sct) <- "RNA"
combined_sct <- NormalizeData(combined_sct, verbose = FALSE)
sct_markers <- FindAllMarkers(combined_sct, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.csv(sct_markers, file.path(out_dir, "sct_markers.csv"), row.names = FALSE)
```

进行差异基因分析（Marker 寻找）。这是一个关键点：虽然我们使用 SCT 进行整合和聚类，但在寻找差异基因时，通常建议回到原始的 counts 数据。

-   `DefaultAssay(combined_sct) <- "RNA"`：将默认 assay 切换回原始 RNA 数据层。
-   `NormalizeData`：因为之前的 SCT 流程跳过了标准的 NormalizeData，所以在进行差异分析前，需要对 "RNA" assay 进行标准的对数标准化（LogNormalize），以确保 FoldChange 计算的解释性符合常规认知。
-   `FindAllMarkers`：寻找每个 cluster 的标记基因。
    -   `only.pos = TRUE`：只保留上调基因。
    -   `logfc.threshold = 0.25`：这里设定的阈值（约 1.19 倍）比 CCA 流程稍低，旨在捕获更多细微差异的基因。

#### Harmony v4

Harmony 是一种基于投影的高效整合算法，它的核心思想并非改变原始的基因表达矩阵，而是在主成分空间（PCA space）中对细胞坐标进行校正。Harmony 通过迭代算法，在最大化不同批次细胞混合（Maximum Diversity Clustering）的同时，最小化同种细胞类型内部的距离，从而去除批次效应。相比 CCA/MNN，Harmony 计算速度更快，内存占用更低，非常适合大规模数据集。

```R
rm(list = ls())

library(Seurat)
library(dplyr)
library(harmony)
library(ggplot2)

options(future.globals.maxSize = 100 * 1024^3)

setwd('/home/jiaoyuan/01.projects/pal-scr')

out_dir <- "results"
if (!dir.exists(out_dir)) {
  dir.create(out_dir)
}

bud <- readRDS("bud/bud.RDS")
cytoledon <- readRDS("cytoledon/cytoledon_hypocotyl.RDS")
leaf <- readRDS("leaf/leaf.RDS")
root <- readRDS("root/root.RDS")
stem1 <- readRDS("stem/pal0w1.RDS")
stem2 <- readRDS("stem/pal0w2.RDS")
stem3 <- readRDS("stem/pal0w3.RDS")

stem <- merge(stem1, y = c(stem2, stem3))

bud$sample <- 'bud'
cytoledon$sample <- 'cytoledon'
leaf$sample <- 'leaf'
root$sample <- 'root'
stem$sample <- 'stem'
```

环境初始化与数据载入。这里首先处理了 `stem` 的生物学重复，因为它们来自同一组织类型，理应不存在强烈的批次效应，直接使用 `merge` 进行物理合并。

```R
bud <- RenameCells(bud, add.cell.id = 'bud')
cytoledon <- RenameCells(cytoledon, add.cell.id = 'cytoledon')
leaf <- RenameCells(leaf, add.cell.id = 'leaf')
root <- RenameCells(root, add.cell.id = 'root')
stem <- RenameCells(stem, add.cell.id = 'stem')

scRNA_harmony <- merge(x = bud, y = c(cytoledon, leaf, root, stem))

rm(bud, cytoledon, leaf, root, stem1, stem2, stem3, stem)
gc()
```

数据合并与清理。与 CCA/SCT 流程一开始保留为列表（list）不同，Harmony 流程通常要求先将所有数据合并为一个巨大的 Seurat 对象。

-   `RenameCells`：给细胞条形码（Barcode）加上样本前缀，防止不同测序批次中出现重复的 Barcode 导致冲突。
-   `merge`：将所有对象合并为 `scRNA_harmony`。
-   `rm` 和 `gc`：删除原有的独立对象并强制进行垃圾回收，释放内存，这对处理大对象至关重要。

```R
scRNA_harmony <- NormalizeData(scRNA_harmony, normalization.method = "LogNormalize", scale.factor = 10000)
scRNA_harmony <- FindVariableFeatures(scRNA_harmony, selection.method = "vst", nfeatures = 2000)
```

对合并后的对象进行标准的预处理。

-   `NormalizeData`：进行全局对数标准化：
    $$
    \ln \left( \frac{count}{TotalCounts} \times 10000 + 1 \right)
    $$
-   `FindVariableFeatures`：寻找 2000 个高变基因。注意这里是在合并后的数据集中寻找高变基因，意味着这些基因主要捕捉的是整体的生物学差异，但也可能包含批次差异（Harmony 稍后会处理这个问题）。

```R
s_genes <- read.table('pal_S.txt', header = FALSE, stringsAsFactors = FALSE)$V1
g2m_genes <- read.table('pal_G2M.txt', header = FALSE, stringsAsFactors = FALSE)$V1

scRNA_harmony <- CellCycleScoring(scRNA_harmony,
                                  s.features = s_genes,
                                  g2m.features = g2m_genes,
                                  set.ident = TRUE)

scRNA_harmony <- ScaleData(scRNA_harmony,
                           vars.to.regress = c("S.Score", "G2M.Score"))
```

细胞周期评分与回归。

-   `CellCycleScoring`：计算 S 期和 G2M 期得分。
-   `ScaleData`：这是关键的一步。在进行 Z-score 标准化的同时，指定 `vars.to.regress = c("S.Score", "G2M.Score")`。这会构建线性模型，将细胞周期得分作为协变量，计算出的残差（Residuals）将作为后续 PCA 分析的输入。这意味着进入 PCA 的数据已经剔除了细胞周期造成的影响。

```R
scRNA_harmony <- RunPCA(scRNA_harmony, npcs = 50, verbose = FALSE)

scRNA_harmony <- RunHarmony(scRNA_harmony, group.by.vars = "sample")
```

运行 Harmony 整合。

-   `RunPCA`：首先必须计算主成分（PCA），因为 Harmony 是在 PC 空间操作的。这里计算了前 50 个 PC。
-   `RunHarmony`：这是核心函数。
    -   `group.by.vars = "sample"`：告诉算法依据哪一列元数据来校正批次效应（这里是样本来源）。
    - 原理：Harmony 并不修改基因表达矩阵（assay 数据不变），而是创建了一组新的、经过校正的“嵌入（Embeddings）”，存储在 `scRNA_harmony@reductions$harmony` 中。它通过软聚类（Soft Clustering）将相似的细胞拉近，利用线性校正将属于同一簇但不同批次的细胞中心对齐。

```R
scRNA_harmony <- FindNeighbors(scRNA_harmony, reduction = "harmony", dims = 1:30)
scRNA_harmony <- FindClusters(scRNA_harmony, resolution = 0.5)

scRNA_harmony <- RunUMAP(scRNA_harmony, reduction = "harmony", dims = 1:30, reduction.name = "umap_harmony")
```

基于整合后的数据进行聚类和降维。

- 关键点：`reduction = "harmony"`。所有的下游分析（构图、聚类、UMAP）都必须显式指定使用 Harmony 校正后的嵌入，而不是原始的 "pca"。
-   `dims = 1:30`：使用前 30 个 Harmony 维度。

```R
plot_umap <- DimPlot(scRNA_harmony, reduction = "umap_harmony", label=T, raster=FALSE)
ggsave(plot=plot_umap, file.path(out_dir, 'harmony_umap.pdf'), width=8, height=7)

plot_umap_sample <- DimPlot(scRNA_harmony, reduction = "umap_harmony", group.by = "sample", raster=FALSE)
ggsave(plot=plot_umap_sample, file.path(out_dir, 'harmony_umap_sample.pdf'), width=8, height=7)
```

可视化结果。

- 生成 `harmony_umap.pdf`：展示聚类结果。
- 生成 `harmony_umap_sample.pdf`：展示样本分布。如果 Harmony 整合成功，不同样本的细胞应该在 UMAP 图上充分混合（即不再按 Sample 分群），除非存在真实的生物学特异性细胞类型。

```R
harmony_markers <- FindAllMarkers(scRNA_harmony, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.csv(harmony_markers, file.path(out_dir, "harmony_markers.csv"))

saveRDS(scRNA_harmony, file = "combined_harmony.rds")
```

寻找差异基因并保存。

-   `FindAllMarkers`：虽然前面的聚类是基于 Harmony 嵌入的，但差异分析是回到默认的 Assay（通常是 "RNA"）中的 `data` 层（即 LogNormalized counts）进行的。这是 Harmony 的优势之一：它保留了原始表达值的物理意义，不像 CCA/MNN 可能会生成修正后的表达矩阵（Integrated Assay），因此直接在原始标准化数据上找 Marker 更加直观和准确。

### UMAP 参数



### 重聚类

### 拟时序分析
