# 组学分析流程合集 (omics-pipelines)

本仓库收集了若干基于 [Snakemake](https://snakemake.readthedocs.io/) 编写、并用
[pixi](https://pixi.sh/) 管理软件环境的二代测序数据处理流程。每个流程都是一个
独立的子目录，输入是原始的双端测序 FASTQ 文件，输出是可直接用于下游分析的结果
（QC 报告、比对 BAM、bigWig 覆盖度轨迹、定量矩阵、差异分析表等）。

| 子目录 | 流程 | 简介 | 详细文档 |
|--------|------|------|----------|
| [`rnaseq/`](rnaseq/README.md)   | RNA-seq    | fastp → STAR 比对 → 基因计数 → DESeq2 差异表达 | [README](rnaseq/README.md) |
| [`chipseq/`](chipseq/README.md) | ChIP-seq   | bowtie2 比对 → 过滤 → MACS3 call peak（含 input 对照）→ 注释 | [README](chipseq/README.md) |
| [`atacseq/`](atacseq/README.md) | ATAC-seq   | bowtie2 比对 → 过滤 + Tn5 位移 → MACS3 call peak → 注释 | [README](atacseq/README.md) |
| [`reseq/`](reseq/README.md)     | 重测序(DNA) | fastp → bwa mem → 标记重复 → 分析就绪 BAM（不含变异检测） | [README](reseq/README.md) |

每个子目录里还有一份**更详细的英文 README**，包含完整的流程图、参数说明和输出文件
清单。本文件只讲解所有流程**通用**的安装与运行方式。

---

## 一、需要自己安装 Snakemake 吗？

**不需要。** 你**只需要安装 pixi 这一个工具**即可。

- Snakemake 本身被声明在每个流程的 `pixi.toml` 里，pixi 会把它装进一个隔离的
  `snakemake` 环境，运行时用 `pixi run -e snakemake snakemake ...` 调用。
- 流程中用到的其它生信软件（STAR、bowtie2、bwa、samtools、MACS3、deepTools、
  R/DESeq2 等）同样由 pixi 自动安装到各自隔离的环境里，**你都不用手动装**。

也就是说，整个工具链里你唯一要手动安装的只有 pixi。

### 安装 pixi

```bash
curl -fsSL https://pixi.sh/install.sh | bash
```

安装完成后重新打开终端（或 `source` 你的 shell 配置），确认可用：

```bash
pixi --version
```

> 运行环境要求：Linux x86-64（所有 `pixi.toml` 都固定了
> `platforms = ["linux-64"]`）。RNA-seq 的 STAR 建索引对内存要求较高（完整人类
> 基因组约需 ~30 GB 内存）。

---

## 二、通用使用流程

下面以 `rnaseq` 为例，其它流程把目录名换掉即可，步骤完全一致。

### 1. 进入流程目录并安装环境（仅第一次）

```bash
cd rnaseq

# 安装该流程声明的所有 pixi 环境（会下载工具，耗时较长，只需一次）
pixi install
```

> ⚠️ ChIP-seq / ATAC-seq 的 ChIPseeker 注释步骤还需要**手动**安装物种注释包，
> 详见各自 README 的 “Setup” 一节（进入 `pixi run -e chipseeker R` 后用
> `BiocManager::install(...)` 安装）。

### 2. 准备参考基因组索引（不由流程生成）

每个流程都需要事先建好的参考索引，索引路径在 `config.yaml` 的 `ref:` 下配置。
不同流程用的比对软件不同，建索引命令也不同：

```bash
# RNA-seq：STAR 索引
pixi run -e star STAR --runMode genomeGenerate \
  --genomeDir ref/star_index --genomeFastaFiles genome.fa \
  --sjdbGTFfile ref/annotation.gtf --sjdbOverhang 100 --runThreadN 16

# ChIP-seq / ATAC-seq：bowtie2 索引
pixi run -e default bowtie2-build genome.fa ref/human

# 重测序：bwa 索引
pixi run -e default bwa index ref/genome.fa
```

### 3. 填写样本表 `samples.tsv`

仓库里自带的 `samples.tsv` 只是**占位示例**，请替换成你自己的样本。它是
**制表符（Tab）分隔**的表格，每行一个样本。各流程的列略有不同：

- **RNA-seq**：`sample / group / layout / r1 / r2`，其中 `group` 决定 DESeq2 的
  比较分组——`samples.tsv` 中**第一个出现的 group 作为参照（WT）**，其余每个 group
  都与它比较。
- **ChIP-seq**：在上面基础上多了 `type`（IP/input）和 `control` 列，`control` 填
  对应 input 样本的 `sample` 名，用来驱动 MACS3 的 `-c` 对照。
- **ATAC-seq / 重测序**：`sample / group / layout / r1 / r2`。

`r1` / `r2` 是双端 FASTQ 的路径（相对于该流程目录）。

### 4. 按需调整参数 `config.yaml`

参考路径、各软件线程数、fastp/比对/call peak 等参数都集中在 `config.yaml` 里，
可按需修改（例如 RNA-seq 的链特异性 `star.strandedness`、ChIP-seq 的窄峰/宽峰等）。
更换物种时也是改这里（详见各 README 的 “Adapting to another organism”）。

### 5. 运行流程

Snakemake 在它自己的 pixi 环境里运行，每条规则再去调用对应工具的环境：

```bash
# 干跑（dry-run）：只打印将要执行的步骤，不真正运行，强烈建议先跑一次
pixi run -e snakemake snakemake -n

# 正式运行整个流程（-j 后是并行核数，按机器调整）
pixi run -e snakemake snakemake -j 32 -p

# 只生成某个目标文件
pixi run -e snakemake snakemake -j 8 results/04-quant/counts.tsv

# 可视化任务依赖图(DAG)
pixi run -e snakemake snakemake --dag | dot -Tsvg > dag.svg
```

结果默认输出到各流程目录下的 `results/`，每条规则的日志写在 `logs/`。

---

## 三、小结

1. 装 **pixi**（唯一需要手动安装的工具）。
2. `cd <流程目录>` → `pixi install`。
3. 建好参考索引，改 `samples.tsv` 和 `config.yaml`。
4. `pixi run -e snakemake snakemake -n` 先干跑，再去掉 `-n` 正式运行。

更多细节请查看每个子目录里的 README。
