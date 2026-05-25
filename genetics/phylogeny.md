# Phylogenetic Analysis

## 安装环境

```bash
mamba create -n pa -c bioconda -c conda-forge paml mafft iqtree figtree gotree biopython python=3.9 -y
mamba activate pa
```

## 构建进化树

常见的进化树算法如下：

- 邻接法（Neighbor-Joining）：基于距离矩阵的聚类算法，通过反复连接总进化支长最小的分类单元构建无根树，计算速度快、适合处理上千条序列的大数据集，在序列分歧较大时易受长枝吸引影响而导致拓扑结构错误，适用于序列相似度较高、进化速率均匀的快速建树或大规模数据的初步建树。
- 最大似然法（Maximum Likelihood）：在预设的氨基酸替代模型下，搜索能使观测数据出现概率最大的树拓扑结构及其枝长，统计基础扎实、准确性高且能进行模型比较和分支支持率评估，计算量极大、对大数据集耗时过长，适用于对准确性要求严格、序列数量适中（通常建议少于 200 条）且分歧度较高的精细系统发育研究。
- 最大简约法（Maximum Parsimony）：基于奥卡姆剃刀原则，寻找能够用最少氨基酸替换次数解释序列差异的树拓扑结构，不依赖复杂的进化模型、避免模型误用，且在序列差异较小的近缘物种分析中表现稳健，当序列分歧大或存在多位点平行突变时易产生长枝吸引和统计不一致，适用于序列分歧度低、位点突变较为简单的近缘物种或小规模同源蛋白分析。
- 贝叶斯推断法（Bayesian Inference）：利用马尔可夫链蒙特卡洛（MCMC）算法对树空间进行采样，根据后验概率分布评估不同树拓扑结构和模型参数的可能性，能直接给出节点支持率的概率解释，可灵活整合先验知识，结果直观可靠，计算资源消耗大、对先验设定敏感且收敛诊断较为复杂，适用于需要高置信度分支支持、序列数量中等（通常少于 100 条）且希望结合先验生物学信息的精细系统发育推断。

先使用蛋白质进行序列比对：

```bash
mafft --auto genes.fasta > genes_aligned.aln
```

使用 Trimal 去除比对质量差的区域（杂乱区），提高建树的准确性：

```bash
trimal -in genes_aligned.aln -out genes_aligned_trimmed.aln -automated1
```

使用 IQ-TREE 进行建树：

```bash
iqtree -s genes_aligned_trimmed.aln -m MFP -bb 1000 -nt AUTO
```

- `-s` 指定序列
- `-m MFP` 自动选择最佳模型
- `-bb 1000` 快速自展法
- `-nt AUTO` 使用多核
- `genes_aligned.fasta.treefile` 是 Newick 格式的进化树

也可以在建树的同时定根：

```bash
iqtree -s genes_aligned_trimmed.aln -m MFP -bb 1000 -o Os03t0122300-01,Os04t0581100-01 -asr -nt AUTO
```

## 重建祖先序列

PAML 建议使用有根树，可以使用 gotree 打开该树文件，手动指定外类群并保存为 `input.tre`：

列出树中所有的物种名（Tips/Leaves）：

```bash
gotree labels -i genes_aligned_trimmed.aln.treefile
```

假设你的外类群名称是 SpeciesA，定根（Reroot）：

```bash
gotree reroot outgroup -i genes_aligned_trimmed.aln.treefile -l "SpeciesA" -o rooted.tree
```

也可以通过指定基因的方式：

```bash
cat > outgroup.txt <<EOF  
AT4G10490.1  
AT4G10500.1  
AT5G24530.1  
EOF
```

- `-i` 输入原始树
- `-l` 指定外类群名字
- `-o` 输出定根后的树

```bash
gotree reroot outgroup -i genes_aligned_trimmed.aln.treefile -l outgroup.txt -o rooted.tree
```

清洗树文件，去除 Bootstrap 等 PAML 不识别的信息，只保留拓扑和支长：

```bash
sed -E 's/\)[0-9\/.]+:/\):/g' rooted.tree > input.tree
```

PAML 对输入格式非常挑剔，通常需要 PHYLIP 格式，可以写一个简单的 Python 脚本（或使用 PhyloSuite 转换）将 FASTA 转换为 PHYLIP。这里使用 PAML 的 coding 模式，也就是使用 CDS 序列来重建祖先序列。

创建一个 `convert.py` 来转换格式

```python
import sys, os, re, subprocess, argparse

def run():
    p = argparse.ArgumentParser()
    p.add_argument("-fa", required=True)
    p.add_argument("-itre", required=True)
    p.add_argument("-phy", required=True)
    p.add_argument("-otre", required=True)
    args = p.parse_args()

    raw, order = {}, []
    with open(args.fa, 'r') as f:
        n, s = "", []
        for l in f:
            if l.startswith(">"):
                if n: raw[n] = "".join(s)
                n = l[1:].strip().split()[0]
                order.append(n)
                s = []
            else: s.append(l.strip())
        if n: raw[n] = "".join(s)

    t = {'ATA':'I','ATC':'I','ATT':'I','ATG':'M','ACA':'T','ACC':'T','ACG':'T','ACT':'T',
         'AAC':'N','AAT':'N','AAA':'K','AAG':'K','AGC':'S','AGT':'S','AGA':'R','AGG':'R',
         'CTA':'L','CTC':'L','CTG':'L','CTT':'L','CCA':'P','CCC':'P','CCG':'P','CCT':'P',
         'CAC':'H','CAT':'H','CAA':'Q','CAG':'Q','CGA':'R','CGC':'R','CGG':'R','CGT':'R',
         'GTA':'V','GTC':'V','GTG':'V','GTT':'V','GCA':'A','GCC':'A','GCG':'A','GCT':'A',
         'GAC':'D','GAT':'D','GAA':'E','GAG':'E','GGA':'G','GGC':'G','GGG':'G','GGT':'G',
         'TCA':'S','TCC':'S','TCG':'S','TCT':'S','TTC':'F','TTT':'F','TTA':'L','TTG':'L',
         'TAC':'Y','TAT':'Y','TAA':'*','TAG':'*','TGA':'*','TGC':'C','TGT':'C','TGG':'W'}

    pep = {}
    for n in order:
        s = raw[n].upper()
        s = s[:len(s)-len(s)%3]
        if s[-3:] in ["TAA","TAG","TGA"]: s = s[:-3]
        raw[n] = s
        pep[n] = "".join(t.get(s[i:i+3], "X") for i in range(0, len(s), 3)).replace("*", "")

    with open("tmp.fa", "w") as f:
        for n in order: f.write(f">{n}\n{pep[n]}\n")
    
    subprocess.run(["mafft", "--auto", "--quiet", "tmp.fa"], stdout=open("tmp.aln", "w"), check=True)

    aln_pep, cur = {}, ""
    with open("tmp.aln", "r") as f:
        for l in f:
            if l.startswith(">"):
                cur = l[1:].strip()
                aln_pep[cur] = ""
            else: aln_pep[cur] += l.strip()

    id_map = {n: (re.sub(r'[^a-zA-Z0-9]', '', n)[:7] + f"_{i+1:02d}")[:10] for i, n in enumerate(order)}

    with open(args.phy, "w") as f:
        f.write(f"{len(order)} {len(list(aln_pep.values())[0])*3}\n")
        for n in order:
            res, pos = [], 0
            for c in aln_pep[n]:
                if c == "-": res.append("---")
                else:
                    res.append(raw[n][pos:pos+3])
                    pos += 3
            seq = re.sub(r'[^ATCGN\-]', 'N', "".join(res).upper())
            f.write(f"{id_map[n].ljust(10)}  {seq}\n")

    with open(args.itre, 'r') as f: tr = f.read()
    for n in sorted(order, key=len, reverse=True): tr = tr.replace(n, id_map[n])
    tr = re.sub(r'\)Node\d+', ')', re.sub(r'\)[0-9]+:', '):', tr))
    with open(args.otre, "w") as f: f.write(tr.strip())

    for tmp in ["tmp.fa", "tmp.aln"]:
        if os.path.exists(tmp): os.remove(tmp)

if __name__ == "__main__":
    run()
```

- 脚本需要安装 `mafft`
- 建树需要高质量的比对区域（去除噪声），但祖先重建通常希望得到全长基因，iqtree 用的树是基于修剪后的，而 convert.py 重新对比全长 CDS 进行 ASR，PAML 会根据树的拓扑结构在全长比对序列上计算概率。

使用方法：

```bash
python convert.py -fa input.cds.fa -itre input.tree -phy paml.phy -otre paml.tre
```


在当前目录下创建一个文本文件，命名为 `codeml.ctl` 参考下面的配置模板：

```text
seqfile = paml.phy
     treefile = paml.tre
      outfile = mlc

        noisy = 9
      verbose = 2
      runmode = 0

      seqtype = 1
    CodonFreq = 2
        model = 0
      NSsites = 0

        icode = 0
    fix_kappa = 0
        kappa = 2
    fix_omega = 0
        omega = 0.1

    RateAncestor = 1

   Small_Diff = .5e-6
    cleandata = 0
```

在 PAML 的 `codeml` 运行中，是否需要额外的模型文件（ `jones.dat`）取决于 `seqtype` 设置：

* 密码子模型 (`seqtype = 1`)：基于 CDS 序列重建祖先序列或计算 $dN/dS$。此时 PAML 使用基于遗传密码子的替代矩阵，不需要外部氨基酸替换矩阵文件。
* 蛋白质模型 (`seqtype = 2`)：基于氨基酸序列重建祖先序列。必须指定一个氨基酸替代模型（如 `model = 2`），并在运行目录下提供相应的矩阵文件（如 `jones.dat`, `wag.dat`, `lg.dat` 等），通常位于 `$CONDA_PREFIX/dat/` 文件夹中。

确保 `paml.phy`, `paml.tre`, `codeml.ctl` 都在同一个文件夹下。

```bash
codeml codeml.ctl
```

运行完成后，你会看到以下文件：
- `mlc`：包含模型估计的参数（如支长、替代速率）
- `rst`：最重要的文件，包含重建的祖先序列
- `rub`：运行过程中的临时数据

在 `rst` 中搜索 `Marginal reconstruction of ancestral sequences`，会看到一句话类似 `Prob distribution at node 26, by site`，说明 node26 是整棵树的最原始祖先或者说根节点。

使用 PhyloSuite 的 "PAML Result Parser" 功能，导出所有节点的祖先序列 FASTA 文件，或者使用一个 Python 脚本：

```python
import re
import os

def process():
    if not os.path.exists("rst"):
        return

    with open("rst", "r") as f:
        lines = f.readlines()

    tree_str = ""
    for i, line in enumerate(lines):
        if "tree with node labels" in line or "TreeView" in line:
            for j in range(i + 1, i + 5):
                if j < len(lines) and "(" in lines[j] and ";" in lines[j]:
                    tree_str = lines[j].strip()
                    break
        if tree_str: break
    
    if not tree_str:
        for line in lines:
            if line.strip().startswith("(") and line.strip().endswith(");") and "#" in line:
                tree_str = line.strip()
                break

    if tree_str:
        tree_str = tree_str.replace("#", "").replace(" ", "")
        with open("paml_itol.tre", "w") as f:
            f.write(tree_str)

    content = "".join(lines)
    section = re.search(r"List of extant and reconstructed sequences(.*?)Overall accuracy", content, re.S)
    if section:
        nodes = re.findall(r"node\s+#(\d+)\s+([ATCG\s\-]+)", section.group(1), re.I)
        with open("paml_ancestors.fa", "w") as f:
            for node_id, seq in nodes:
                clean_seq = "".join(seq.split())
                f.write(f">{node_id}\n{clean_seq}\n")

if __name__ == "__main__":
    process()
```

生成的祖先序列与树文件中的节点相对应，可以通过 iTOL 等工具打开树文件查看节点。