# Collinearity Analysis

- 直系同源（Orthologs）是指由于物种形成（Speciation）事件而产生的基因分化，这些基因源自共同祖先，随着物种的分裂而分别进入不同的物种基因组中。在进化过程中，直系同源基因通常在不同物种间保留了相似或相同的生物学功能，因此它们是跨物种功能推断、基因组比较以及构建物种进化树最可靠的依据，其分化点在时间上严格对应物种的分化时刻。
- 旁系同源（Paralogs）则是由于同一物种内部发生的基因复制（Gene Duplication）事件而产生的基因拷贝，这些拷贝通常并存于同一个基因组内（或随进化分配到亲缘物种中）。旁系同源基因在产生后，由于遗传冗余，不同的拷贝往往会发生功能分化：有的保留原功能，有的发生功能弱化或丢失，有的则进化出全新的生物学功能（新功能化）。它们是研究基因家族扩张、收缩以及生物适应性进化的核心研究对象。
- 共线性分析（Synteny Analysis）是指通过比较不同物种或同一物种不同染色体区域间基因的排列顺序和物理位置，来揭示基因组进化联系的研究方法。简单来说，如果两个物种的某段染色体上，不仅基因序列相似，而且这些基因的先后排列顺序也高度一致，我们就称这两个区域具有“共线性区块”。

## 共线性分析

在 mamba 中安装 jcvi 及相关工具：

```bash
mamba create -n jcvi -c bioconda -c conda-forge jcvi last bedtools seqkit emboss -y
mamba activate jcvi
```

将 GFF3 转换为简化的 BED 格式，并提取最长转录本以减少冗余：

```bash
# 以拟南芥 (ath) 和水稻 (osa) 为例
python -m jcvi.formats.gff bed --type=mRNA --key=transcript_id ath.gff3 -o ath.bed
python -m jcvi.formats.gff bed --type=mRNA --key=transcript_id osa.gff3 -o osa.bed

# 只保留每个基因最长的转录本
python -m jcvi.formats.bed uniq ath.bed
python -m jcvi.formats.bed uniq osa.bed
```

- 根据实际分析时需要的基因标识符修改 `--key` 参数

利用去重后的 BED 文件从原始 FASTA 中提取对应的 CDS 序列：

```bash
seqkit grep -f <(cut -f 4 ath.uniq.bed) ath.cds.fa | seqkit seq -i > ath.cds
seqkit grep -f <(cut -f 4 osa.uniq.bed) osa.cds.fa | seqkit seq -i > osa.cds
```

确保 BED 文件和 CDS 文件前缀一致。运行 ortholog 命令，此过程会调用 LAST 进行比对并寻找共线性块（Anchors）：

```bash
cp ath.uniq.bed ath.bed
cp osa.uniq.bed osa.bed
python -m jcvi.compara.catalog ortholog ath osa --cscore=.99 --no_strip_names
```

生成的 `.anchors` 文件包含了所有共线性点。为了绘图美观，需要进行筛选：

```bash
# 筛选包含至少30个基因的区块
python -m jcvi.compara.synteny screen --minspan=30 --simple ath.osa.anchors ath.osa.anchors.new
```

### Macrosynteny

创建 `seqids` 文件，定义要在图中显示的染色体编号：

```text
1,2,3,4,5
Chr1,Chr2,Chr3,Chr4,Chr5,Chr6,Chr7,Chr8,Chr9,Chr10,Chr11,Chr12
```

创建 `layout` 文件，定义画布中物种的位置、颜色及连接关系：

```text
# x, y, rotation, ha, va, color, ratio, label
0.5, 0.6, 0, left, center, #7fc97f, 1, Arabidopsis
0.5, 0.4, 0, left, center, #beaed4, 1, Rice
e, 0, 1
```

运行绘图模块生成 `karyotype.pdf` 文件：

```bash
python -m jcvi.graphics.karyotype seqids layout --shadestyle=line --chrstyle=roundrect
```

可以添加参数：

- `--shadestyle=line`：将填充带改为线条连接
- `--chrstyle=roundrect`：将染色体绘制为圆角矩形

### Microsynteny

使用 `mcscan` 命令提取局部区块，`--iter=1` 表示每个参考区域只取一个最佳匹配块：

```bash
python -m jcvi.compara.synteny mcscan ath.bed ath.osa.lifted.anchors --iter=1 -o ath.osa.block
```

编写 local.layout 配置文件，并合并 BED 文件：

```bash
# 合并 BED 文件
python -m jcvi.formats.bed merge ath.bed osa.bed -o local.bed

# 绘图
python -m jcvi.graphics.synteny ath.osa.block local.bed local.layout --genelabelsize 6
```

- `--genelabelsize 6` 显示基因名称，并设置字体大小 
- `--genelabels=ID1,ID2` 只显示感兴趣的基因 ID
- `--glyphstyle=arrow` 将基因绘制成方向箭头

### Microsynteny getting fancy

将 GFF3 转换为简化的 BED 格式，并提取最长转录本以减少冗余并且从原始 FASTA 中提取对应的 CDS 序列：

```bash
python -m jcvi.formats.gff bed --type=mRNA --key=transcript_id ath.gff3 -o ath.bed
python -m jcvi.formats.bed uniq ath.bed
seqkit grep -f <(cut -f 4 ath.uniq.bed) ath.cds.fa | seqkit seq -i > ath.cds

python -m jcvi.formats.gff bed --type=mRNA --key=transcript_id osa.gff3 -o osa.bed
python -m jcvi.formats.bed uniq osa.bed
seqkit grep -f <(cut -f 4 osa.uniq.bed) osa.cds.fa | seqkit seq -i > osa.cds

python -m jcvi.formats.gff bed --type=mRNA --key=transcript_id zma.gff3 -o zma.bed
python -m jcvi.formats.bed uniq zma.bed
seqkit grep -f <(cut -f 4 zma.uniq.bed) zma.cds.fa | seqkit seq -i > zma.cds
```

执行共线性分析，找出物种间成对出现的直系同源基因（Anchors）：

```bash
cp ath.uniq.bed ath.bed && cp osa.uniq.bed osa.bed && cp zma.uniq.bed zma.bed
python -m jcvi.compara.catalog ortholog ath osa --cscore=.99 --no_strip_names
python -m jcvi.compara.catalog ortholog osa zma --cscore=.99 --no_strip_names
```

筛选去除零星、不连续的噪音点，只保留具有显著共线性的区块：

```bash
python -m jcvi.graphics.dotplot ath.osa.anchors
python -m jcvi.graphics.dotplot osa.zma.anchors
python -m jcvi.compara.synteny screen --minspan=30 --simple ath.osa.anchors ath.osa.anchors.simple
```

准备 `seqids` 文件：

```txt
1,2,3,4,5
Chr1,Chr2,Chr3,Chr4,Chr5,Chr6,Chr7,Chr8,Chr9,Chr10,Chr11,Chr12
1,2,3,4,5,6,7,8,9,10
```

准备 `layout` 文件：

```txt
# y, xstart, xend, rotation, color, label, va, bed
0.7, 0.1, 0.8, 0, #7fc97f, Arabidopsis, top, ath.bed
0.5, 0.1, 0.8, 0, #beaed4, Rice, top, osa.bed
0.3, 0.1, 0.8, 0, #fdc086, Maize, top, zma.bed
# edges
e, 0, 1, ath.osa.anchors.simple
e, 1, 2, osa.zma.anchors.simple
```

绘制核型图：

```bash
python -m jcvi.graphics.karyotype seqids layout --chrstyle=roundrect
```

利用 `mcscan` 提取两两物种间的共线性块，并以水稻 (osa) 为桥梁进行 Join：

```bash
python -m jcvi.compara.synteny mcscan ath.bed ath.osa.lifted.anchors --iter=1 -o ath.osa.block
python -m jcvi.compara.synteny mcscan osa.bed osa.zma.lifted.anchors --iter=1 -o osa.zma.block

python -m jcvi.formats.base join ath.osa.block osa.zma.block --noheader | awk '{print $1"\t"$2"\t"$4}' > all.block
```

准备 `all.bed` 文件：

```bash
python -m jcvi.formats.bed merge ath.bed osa.bed zma.bed -o local.bed
```

准备 `all.layout` 文件：

```
# x, y, rotation, ha, va, color, ratio, label
0.5, 0.8, 0, center, top, , 1, Arabidopsis
0.5, 0.5, 0, center, center, , 1, Rice
0.5, 0.2, 0, center, bottom, , 1, Maize
# edges
e, 0, 1
e, 1, 2
```

绘制局部共线性：

```bash
python -m jcvi.graphics.synteny all.block all.bed all.layout --genelabelsize 6 --glyphstyle=arrow
```