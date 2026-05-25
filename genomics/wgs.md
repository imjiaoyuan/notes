# WGS

## VCF 数据格式

查看列标题

```bash
zcat Allchr.rename.snp.vcf.gz | grep "^#CHROM"
```

```txt
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	103C	105A	108C	11-1A	111A	112-1A	113-1A	114-3A	115-4B	116C	117-2C	118-2A	119-1A	12-1A	122-2A	125-1C	126A  127-1A	132-2C	133-1A	137-1A	138A	13A	141A	148-1A	150-2A	153-1A	155-1A	156-1A	158-3C	161-1A	162-1A	164-2A	167A	169-4A	16A	17-5A	170-2C	172-1A	175A	177-1A	178A	179-2A181A	182-2A	183-1A	185A	194-5A	198-1A	199-1A	2-2A	201-1A	205-1A	206-3A	207-3B	208-1A	210-2B	214-2A	217-2A	219-3A	220A	221-1A	222-3A	223-1A	227A	230-1A	233-3A	235-2C	236-3A237-2A	238-2A	23A	24-1A	242-1A	244-3A	246A	250-2A	252-1A	253A	258-4A	260-3A	261-2A	263-3A	265-2A	266-3A	268-5A	269-1A	26A	27-2A	273-3C	275-2A	280-1A	282-3A	283-2A	288-1A290-3A	293-2B	295-2A	296-2A	298A	3-1A	301-1A	302-4C	303-1A	306-1A	307C	308-1A	310-2A	313-1A	314-4A	317-1A	32-1A	321-3B	323-2A	325-1A	326-1A	33-1A	330-2C	331-3A	332A	334-3A336A	338-1A	339-3A	34-1A	340-1A	341-2A	343-1A	344-4A	348-4A	351A	353A	357-1A	358-2A	360-2A	36161A	367-1A	369-2A	36A	378-1-1A	379-2A	380-1A	381-2A	384-1A	386-1A	39-1A 390A	392A	393-1A	394-1A	396A	399-1A	400-3A	402-2A	41-1A	411B	412-1A	413-2A	415-2A	418-2A	421-1A	428-2A	437-1A	438A	440-1A	442B	443-1A	448-4A	44947A	449A	456-1A	469A  473A	476-2A	477-2A	478-1A	482A	486-1A	487A	49-1A	493-2A	494-1A	495-2A	496-30A	497-4A	499A	501A	505A	506-2A	507A	508-1A	511-2A	515-1A	521-3A	525-1A	526-3A	529-5A	53-4A 532A	535A	536A	538A	543-2A	544A	546A	548-1A	54A	550-2A	553-6A	559-1A	560A	563-2A	579-2A	580A	582-2A	584-1A	585-1A	586-1A	59-1A	60A	62-1A	65-1A	66-1A	68C   7-1A	70-1A	74A	76-1A	78-1A	8-2A	80-1A	83-1A	84-1A	85-4A	87A	891-2A	89A	95-1A	97-1A	98-4A	99-6A
```

-   `#CHROM` 染色体名称
-   POS 在染色体上的位置
-   ID 变异位点的 ID（通常是 dbSNP ID）
-   REF 参考基因组碱基
-   ALT 变异碱基
-   QUAL 质量值
-   FILTER 过滤状态
-   INFO 附加信息
-   FORMAT 样本基因型数据的格式说明
- 从 `103C` 开始到最后 `99-6A` 都是样本名称

查看元信息

```bash
zcat Allchr.rename.snp.vcf.gz | grep "^##" | head -20
```

```txt
##fileformat=VCFv4.2
##FILTER=<ID=PASS,Description="All filters passed">
##ALT=<ID=NON_REF,Description="Represents any possible alternative allele not already represented at this location by REF and ALT">
##FILTER=<ID=FS60,Description="FS > 60.0">
##FILTER=<ID=LowQual,Description="Low quality">
##FILTER=<ID=MQ40,Description="MQ < 40.0">
##FILTER=<ID=MQRankSum-12.5,Description="MQRankSum < -12.5">
##FILTER=<ID=QD2,Description="QD < 2.0">
##FILTER=<ID=QUAL30,Description="QUAL < 30.0">
##FILTER=<ID=ReadPosRankSum-8,Description="ReadPosRankSum < -8.0">
##FILTER=<ID=SOR3,Description="SOR > 3.0">
##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Allelic depths for the ref and alt alleles in the order listed">
##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Approximate read depth (reads with MQ=255 or with bad mates are filtered)">
##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Genotype Quality">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=MIN_DP,Number=1,Type=Integer,Description="Minimum DP observed within the GVCF block">
##FORMAT=<ID=PGT,Number=1,Type=String,Description="Physical phasing haplotype information, describing how the alternate alleles are phased in relation to one another; will always be heterozygous and is not intended to describe called alleles">
##FORMAT=<ID=PID,Number=1,Type=String,Description="Physical phasing ID information, where each unique ID within a given sample (but not across samples) connects records within a phasing group">
##FORMAT=<ID=PL,Number=G,Type=Integer,Description="Normalized, Phred-scaled likelihoods for genotypes as defined in the VCF specification">
##FORMAT=<ID=PS,Number=1,Type=Integer,Description="Phasing set (typically the position of the first variant in the set)">
```

-   `##fileformat=VCFv4.2`：使用的是 VCF 版本 4.2 格式
    过滤条件
-   FS60：FS > 60.0 - 链路间性过滤，>60 项无显著链路间
-   NQ40：MQ < 40.0 - 映射质量低
-   QD1：OD < 2.0 - 质量深度比低
-   QUAL39：QUAL < 30.0 - 总体质量偏低
-   SOR3：SOR > 3.0 - 链化值过滤
-   NQRanKSun-12.5：MQRanKSun < -12.5 - 映射质量极和检验异常
-   ReasProsRanKSun-8：ReasProsRanKSun < -8.0 - 读取位置极和检验异常

格式字段

-   GT：基因型 - 如 O/O, O/1, 1/1 表示基因型
-   AD：等位基因深度 - REF 和 ALT 的测序深度
-   DP：总深度 - 该位点的总测序深度
-   GQ：基因型质量 - 基因型呼叫的质量值
-   PL：基因型似然值 - 三种基因型的 Phrend 尺度似然值
-   PGT：物理分型单倍率 - 物理分相信息
-   PTD：物理分相 ID - 分相组标识
-   PS：分相量 - 分相组起始位置

## GWAS

GWAS（全基因组关联分析）是一种通过在全基因组范围内搜寻特定 DNA 变异（通常是单核苷酸多态性，即 SNP）与目标性状（如身高、产量等）之间统计学关联的方法。它能高效地定位控制复杂性状的关键基因或基因组区域，揭示性状背后的遗传架构（即哪些基因在起作用、影响力有多大）。在植物研究中，GWAS 被广泛用于解析产量、抗逆性及生长发育的分子机制。以杨树株高为例，收集成百上千株遗传背景不同的杨树，精确测量每棵树的株高（表型数据），同时通过高通量测序获取它们的全基因组 SNP 变异（基因型数据）；随后利用统计学模型，逐一比对成千上万个 SNP 位点与株高高低的关系。如果发现某些 SNP 位点在长得高的杨树中出现的频率显著高于矮树，就意味着这些位点附近可能存在影响杨树生长的核心基因，从而为培育速生高产的林木品种提供精准的分子标记。

GWAS 与全基因组选择育种是发现与应用的关系，互为补充。GWAS 旨在“找靶点”，即解析控制株高的关键基因位点，揭示性状形成的遗传原理；而全基因组选择（GS）则侧重于“算总分”，即综合利用全基因组所有变异信息来预测树苗的生长潜力，实现直接选优。二者相辅相成：GWAS 鉴定出的重要位点可作为关键信息输入 GS 模型，从而显著提升对育种材料未来表现的预测准确性。通过这种策略，先由 GWAS 发掘影响长高的遗传基础，再由 GS 高效筛选出具有最优基因组合的个体，能够极大地加速杨树良种的选育进程。

目前在实际应用中最常用、最核心的 SNP 类型是 连锁不平衡中的标签 SNP。

### Pipeline

在分析之前，首先要对表型进行标准化。使用 load_traits 函数处理表型性状。表型数据需检查正态性，对偏态分布进行适当转换，确保符合线性模型假设。

样本匹配阶段使用 `bcftools query -l` 获取 VCF 文件中的样本名，并与表型样本名取交集。生成的 `keep.txt` 包含保留样本的原始 ID，`rename.txt` 记录原始 ID 到清洗后 ID 的映射关系。随后的 `bcftools` 处理步骤通过管道符串联执行：

```python
subprocess.run("bcftools view -S keep.txt VCF_FILE | bcftools reheader -s rename.txt | bcftools view -O z -o clean.vcf.gz", shell=True, check=True)
```

-   `-S`：指定需要保留的样本名单列表文件。
-   `-s`：指定包含“旧 ID 新 ID”对应关系的映射文件，用于重命名 VCF 表头。
-   `-O z`：设置输出格式为压缩的 VCF（.vcf.gz）。

`clean.vcf.gz` 是用于后续分析的清洗后基因型文件。

然后对基因型数据进行一下数据过滤，减少假阳性结果。需要过滤低质量 SNP 和个体，包括缺失率过高、稀有变异和偏离哈迪 - 温伯格平衡的位点：

```python
subprocess.run("plink --vcf clean.vcf.gz --geno 0.1 --mind 0.1 --maf 0.05 --hwe 1e-6 --make-bed --out genotype_qc --allow-extra-chr --const-fid", shell=True)
```

-   `--geno 0.1`：剔除缺失率高于 10% 的 SNP 位点
-   `--mind 0.1`：剔除基因型缺失率高于 10% 的个体
-   `--maf 0.05`：过滤次要等位基因频率低于 5% 的稀有变异
-   `--hwe 1e-6`：剔除严重偏离哈迪 - 温伯格平衡的位点（p<1e-6）

基因型格式转换利用 `PLINK` 工具将文本 VCF 转换为二进制格式，提高关联分析的计算效率：

```python
subprocess.run("plink --vcf clean.vcf.gz --make-bed --out genotype --allow-extra-chr --const-fid --silent", shell=True)
```

-   `--vcf`：输入待转换的 VCF 基因型文件。
-   `--make-bed`：输出 PLINK 二进制格式文件（.bed, .bim, .fam）。
-   `--allow-extra-chr`：允许处理杨树基因组中非标准的染色体编号（如 Scaffold 名称）。
-   `--const-fid`：强制将所有个体的家系 ID（Family ID）设为 0，适用于植物自然群体。

随后通过 `PLINK` 执行群体结构主成分分析（PCA），提取群体的演化背景信息作为协变量：

```python
subprocess.run("plink --bfile genotype --pca 3 --out pca", shell=True)
```

-   `--bfile`：读取 PLINK 格式的二进制基因型文件前缀。
-   `--pca 3`：计算全基因组前 3 个主成分，用于反映群体分层。
-   `--out`：指定输出文件的前缀。

亲缘关系矩阵由 `GEMMA` 软件计算，其数学公式为：$$K = \frac{ZZ'}{2\sum p_j(1-p_j)}$$ 该矩阵通过捕捉减数分裂中染色体重组产生的遗传相似性，用于校正模型中的随机效应：

```python
subprocess.run("gemma -bfile genotype -gk 1 -o kinship -outdir output", shell=True)
```

-   `-gk 1`：计算中心化相关的亲缘关系矩阵（Centered Relatedness Matrix）。
-   `-outdir`：指定计算结果的存储目录。

最后执行线性混合模型（LMM）关联分析扫描：$$y = Xb + Zg + e$$

```python
subprocess.run("gemma -bfile genotype -k output/kinship.cXX.txt -lmm 1 -p trait.pheno -c covariates.txt -o trait_name -outdir output", shell=True, check=True)
```

-   `-lmm 1`：指定使用 Wald 检验算法进行关联分析。
-   `-k`：读入预先计算好的亲缘关系矩阵文件（Kinship）。
-   `-p`：读入包含杨树株高或叶形 PC 值的表型文件。
-   `-c`：读入包含 PCA 特征向量的协变量文件，用于校正固定效应（群体分层）。
-   `-o`：指定输出结果文件的名称。

### 可视化

- 曼哈顿图：用于展示全基因组关联分析（GWAS）中各个位点与性状的关联显著性。横坐标为染色体位置，纵坐标为‑log₁₀(P 值)，点越高表示关联越强。图中常设置一条显著性阈值线（如 P<5×10⁻⁸），超过阈值的位点被视为潜在关联信号。若为数量性状，成簇的显著点比单一点更可靠，可降低假阳性风险。
-   QQ 图：用于评估 GWAS 结果中 P 值分布的合理性。将实际观测的 P 值与理论期望的 P 值进行比较，点大多沿对角线分布说明结果无严重偏差；若尾部向上偏离对角线，表明存在真实关联信号；若整体向下偏离，则提示可能存在群体结构等系统性误差。该图是判断假阳性和假阴性的重要工具。
-   LD-Block 图：用于显示单核苷酸多态性（SNP）之间的连锁不平衡（LD）关系。通常以三角热图形式呈现，颜色深浅代表 LD 强度（如深色表示强连锁）。该图可帮助识别基因组中共同遗传的 SNP 区域，辅助定位功能变异，并解释多个相邻 SNP 对性状的协同影响。
