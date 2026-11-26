# Marine：TaxDistill物种分类校正

环境：`taxnometer_ntnf_KD`。参数：**α=0.5，T=4，teacher=500M，seed=24，batch=64，epochs=100**。

## 运行命令

以Marine的MMseqs2初始分类结果进行校正：

```bash
conda activate taxnometer_ntnf_KD
cd /disk1/yerongye/TaxDistill/KD_Taxnometer_ntnf/vamb

DATA_DIR=/disk1/yerongye/TaxDistill/Taxnometer/vamb/data/Marine
RUN_DIR=/disk1/yerongye/TaxDistill/KD_Taxnometer_ntnf/vamb/output_marine
mkdir -p "$RUN_DIR"

CUDA_VISIBLE_DEVICES=3 PYTHONPATH="$PWD" python vamb/__main__.py taxometer \
  --outdir "$RUN_DIR/mmseq2_alpha05_T4_seed24" \
  --fasta "$DATA_DIR/contigs_marine.fna.gz" \
  --abundance_tsv "$DATA_DIR/abundances.tsv" \
  --taxonomy "$DATA_DIR/mmseq2_marine_taxonomy.tsv" \
  --genomeocean /disk1/yerongye/TaxDistill/Taxnometer/vamb/500M \
  --kd_alpha 0.5 --kd_temp 4.0 -pt 64 -pe 100 --seed 24 --cuda
```

校正Kraken2或Metabuli结果时，更换`--taxonomy`并使用新的输出目录：

| 初始分类器 | Marine标签文件（位于上述`DATA_DIR`） |
|---|---|
| MMseqs2 | `mmseq2_marine_taxonomy.tsv` |
| Kraken2 | `kraken2_Marine_output_gtdb.tsv` |
| Metabuli | `metabuli_results_taxometer.tsv` |

## Kraken2标签转换后再输入TaxDistill

流程：**Kraken2原始`--output`结果 → NCBI taxid映射到GTDB R226 → 两列标签TSV → TaxDistill的`--taxonomy`**。

- 原始Marine结果：`/disk1/yerongye/TaxDistill/Taxnometer/vamb/data/Marine/kraken2_Marine_output.txt`。
- 转换脚本：`/disk1/yerongye/TaxDistill/Taxnometer/vamb/eval/kraken2_to_gtdb.py`。
- 映射表：`/disk1/yerongye/TaxDistill/Taxnometer/vamb/eval/bac120_metadata_r226.tsv`（这个可以直接让GPT下载），按`ncbi_taxid`查找`gtdb_taxonomy`。

## 核心输入文件：是什么、怎么得到

| 输入参数 | 文件 | 内容 | 获取方式 |
|---|---|---|---|
| `--fasta` | `contigs_marine.fna.gz` | 组装后的contig序列 | 原始reads经质控、组装得到contigs；多样本组装结果可用`Taxnometer/vamb/src/concatenate.py`合并并统一ID。本次直接使用已有Marine文件。 |
| `--abundance_tsv` | `abundances.tsv` | 每条contig在各样本中的丰度；首列`contigname`，后面每个样本一列，本Marine文件有10个样本列 | 各样本reads使用`strobealign --aemb`分别回比对到同一份合并contig FASTA，再用`Taxnometer/vamb/src/merge_aemb.py`合并。 |
| `--taxonomy` | 上表三个TSV之一 | 初始分类标签，表头为`contigs`、`predictions`，谱系标签含物种及上级分类 | 用MMseqs2、Kraken2或Metabuli分类同一份contig FASTA，再转换为两列TSV。Kraken2原始NCBI标签需先映射到GTDB R226。 |

三类数据的contig ID必须对应；标签文件须覆盖参与校正的contig，未分类项可为空。原始reads用于组装和计算丰度，不直接传入TaxDistill。

`--genomeocean`另需已有的GenomeOcean预训练模型目录：`/disk1/yerongye/TaxDistill/Taxnometer/vamb/500M`，其中包含权重、`config.json`和tokenizer文件；使用发布的预训练checkpoint或复制已有完整目录。

## 三个分类器：当前数据库、版本和分类参数

| 工具 | 软件版本 | 数据库 | 主要分类参数 |
|---|---|---|---|
| Metabuli | 1.1.1 | GTDB R226 | `--seq-mode 3 --threads 32 --max-ram 256 --lineage 1` |
| MMseqs2 | 18.8cc5c | GTDB R226蛋白库 | `--tax-lineage 1 --threads 32 --split-memory-limit 256G --remove-tmp-files 1` |
| Kraken2 | 2.17.1 | Kraken2_small（NCBI） | `--threads 32 --minimum-hit-groups 3 --confidence 0` |

数据库位置：

```text
Metabuli：/disk1/yerongye/TaxDistill/Taxnometer/vamb/data/databases/Metabuli_GTDB
MMseqs2： /disk1/yerongye/TaxDistill/Taxnometer/vamb/data/databases/GTDB
Kraken2： /disk1/yerongye/TaxDistill/Taxnometer/vamb/data/Gastrointenstinal/Kraken2_small/databases/kraken2_db
```

## 三个数据库的下载命令

以下只记录下载方法；**Metabuli和MMseqs2下载的是固定R226的建库数据，下载后仍需建立各自索引**。

```bash
unset http_proxy https_proxy all_proxy HTTP_PROXY HTTPS_PROXY ALL_PROXY ftp_proxy FTP_PROXY
DB_DOWNLOAD_ROOT=/disk1/yerongye/TaxDistill/database_downloads
R226_BASE=https://data.gtdb.ecogenomic.org/releases/release226/226.0
ARIA2C=(/disk1/yerongye/miniforge3/envs/blast/bin/aria2c --no-conf=true \
  --all-proxy= --http-proxy= --https-proxy= --ftp-proxy= -x 8 -s 8 -k 1M -c)
```

### Metabuli：GTDB R226核酸参考

当前本机库由完整R226代表基因组自建；以下是相同版本的核酸原料下载。Metabuli 1.1.1内置`databases GTDB`列出的库为R214.1，不能用该命令当作R226下载。[GTDB官方R226目录](https://data.gtdb.ecogenomic.org/releases/release226/226.0/genomic_files_reps/)

```bash
"${ARIA2C[@]}" --dir="$DB_DOWNLOAD_ROOT/metabuli_R226" \
  --checksum=md5=b4e6e844bb1f44bcaa13fd88ef1d67d0 \
  "$R226_BASE/genomic_files_reps/gtdb_genomes_reps_r226.tar.gz"
```

建库还需本机已有的`GTDB_TAXDUMP/gtdb-taxdump-R226/`及其中`taxid.map`；建库参数为`--gtdb 1 --threads 64 --max-ram 128`。

### MMseqs2：GTDB R226蛋白参考与分类谱系

固定下载R226，避免`mmseqs databases GTDB`获取随时间变化的`latest`。[MMseqs2 18工作流源码](https://raw.githubusercontent.com/soedinglab/MMseqs2/18-8cc5c/data/workflow/databases.sh)

```bash
"${ARIA2C[@]}" --dir="$DB_DOWNLOAD_ROOT/mmseqs2_R226" \
  --checksum=md5=07577d0438fdee8a64b385b7c38d4c24 \
  "$R226_BASE/genomic_files_reps/gtdb_proteins_aa_reps_r226.tar.gz"

"${ARIA2C[@]}" --dir="$DB_DOWNLOAD_ROOT/mmseqs2_R226" \
  "$R226_BASE/bac120_taxonomy_r226.tsv" \
  "$R226_BASE/ar53_taxonomy_r226.tsv" \
  "$R226_BASE/VERSION.txt"
```

上述GTDB压缩包校验值来自[官方MD5清单](https://data.gtdb.ecogenomic.org/releases/release226/226.0/MD5SUM.txt)。

### Kraken2：公开小型预建库下载示例(可以换成大库)

```bash
"${ARIA2C[@]}" --dir="$DB_DOWNLOAD_ROOT/kraken2_standard16_20261022" \
  --checksum=md5=98dfe7527064392b05b17999c9c084bc \
  https://genome-idx.s3.amazonaws.com/kraken/k2_standard_16_GB_20261022.tar.gz

"${ARIA2C[@]}" --dir="$DB_DOWNLOAD_ROOT/kraken2_standard16_20261022" \
  https://genome-idx.s3.amazonaws.com/kraken/standard_16_GB_20261022/standard_16_GB.md5

cd "$DB_DOWNLOAD_ROOT/kraken2_standard16_20261022"
tar -xzf k2_standard_16_GB_20261022.tar.gz
md5sum -c standard_16_GB.md5
```

## PyTorch与Hugging Face包版本

| 包 | 版本 |
|---|---|
| PyTorch | 2.12.0.dev20260306+cu128 |
| transformers | 4.39.3 |
| huggingface-hub | 0.36.2 |
| tokenizers | 0.15.2 |
| safetensors | 0.7.0 |
| peft | 0.10.0 |
| accelerate | 1.13.0 |
