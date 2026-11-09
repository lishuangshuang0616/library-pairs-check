# library-pair-check

用于 **DNBelab C4 多样本 cDNA / oligo 原始 R1 FASTQ** 的 barcode 指纹配对检查，Rust v0.3 目录版。

输入沿用 `dnbc4tools_py` 的 `--fastqs` 目录约定：每个样本一个目录，里面有 `cDNA/` 与 `oligo/` 两个 FASTQ 子目录。工具自动发现所有样本，把每个样本的 cDNA 与**同一样本及其它所有样本**的 oligo 交叉比较，检查目录约定的配对是否得到 barcode 证据支持，并给出疑似跨样本错配 / 样本交换的线索。**不能据此证明同一供者、同一生物学样本或同一批细胞，也没有经过真实数据校准的自动匹配阈值。**

与 `dnbc4tools_py` 的关系：目录别名、FASTQ 扩展名和 R1/R2 文件名识别按它的 `infrastructure/filesystem/fastq_detector.py`、`fastq_naming.py` 约定实现；它一次只处理一个 `--fastqs` 根、样本名来自 `--name`，且 cDNA↔oligo 的配对只是目录约定，没有 barcode 层面的校验。本工具把多个样本根并列扫描、用目录名作样本名，并补上缺失的文库配对判定。

## 构建

```bash
cargo build --release --locked
./target/release/library-pair-check --help
```

只需要 R1；不读取 R2 内容、不运行 STAR，也不做基因表达、UMI 分子交集或磁珠合并关系验证。支持标准四行 FASTQ、CRLF、普通文本和 gzip（根据文件头识别，支持 concatenated gzip）。不支持换行包装的 FASTQ。

## 输入目录

```
/data/project/                       # --input 指向这里
├── SampleA/
│   ├── cDNA/                        # 也接受 cdna / CDNA
│   │   ├── SampleA_L001_R1_001.fastq.gz
│   │   ├── SampleA_L001_R2_001.fastq.gz
│   │   └── SampleA_L002_R1_001.fastq.gz
│   └── oligo/                       # 也接受 Oligo / OLIGO
│       ├── SampleA_L001_R1_001.fastq.gz
│       └── SampleA_L001_R2_001.fastq.gz
├── SampleB/
│   ├── cDNA/...
│   └── oligo/...
└── ...
```

- `--input DIR` 可重复：每个值可以是「含多个样本目录的根」，也可以直接是「单个样本目录」；样本名取目录名，必须唯一。
- 每个样本目录必须同时存在 `cDNA/` 与 `oligo/`；子目录只扫描一层，文件按名称排序。
- FASTQ 扩展名：`.fastq.gz`、`.fq.gz`、`.fastq`、`.fq`。R1 / R2 由**目录项文件名**的 token 判定（`_R1_`、`_R1`、`_1`、`R1`、`_read1` 及对应 R2），同时命中时取位置靠后的 token；单个无法判定的 FASTQ 视为 R1，多个无法判定则报错。
- 允许符号链接（集群数据常见）：按链接名判断 R1/R2，报告与去重使用解析后的规范路径，链接目标文件名可以不同；悬空链接跳过。若两个文库的链接指向同一真实文件，会因规范路径重复而报错。
- 同一文库目录下所有 R1 文件（多个 lane）合并为一个文库 profile；R2 文件只记录在报告中、不读取，也不校验 R1/R2 是否成对。
- 同一个规范路径不能在多个样本 / 文库中重复出现；硬链接造成的内容重复不检测。输出文件不能覆盖任何输入 FASTQ。

## 检查一组文库

```bash
./target/release/library-pair-check \
  --input /data/project \
  --output pairing-report.json
```

所有 cDNA × 所有 oligo 都会比较，不做强制一对一分配。默认用 **`--sampling reservoir`** 扫描每个文件的全部记录，并以可复现的蓄水池抽样保留最多 100 万条（`--max-records-per-file`）。这会读取完整 gzip 流，不能当作"仅读取前 100 万条"的快速模式；支持 gzip 完整性及文件尾部解析检查。取样包含有效和无效记录，随后独立过滤 barcode，避免只抽高质量记录改变样本分母。

- `--max-records-per-file 0`：保留全部记录，不抽样。证据不足时建议增大上限或完整统计。
- `--sampling prefix`：只读每个文件前 N 条，速度可能更快，但对记录顺序敏感，未读尾部不做完整性检查。
- `--seed 42`：基本随机种子；每个文库、文件使用独立随机流。每文件实际种子写入 JSON，随机流由基本种子、文库角色及规范文件路径确定；增加样本或调换参数顺序不会改变已有文件的抽样。重命名/移动输入路径、修改种子或重排记录会改变具体抽样，不承诺这些情况下结果逐位不变。

每个文件独立检测结构，并独立执行记录上限；同一文件内部应保持结构一致。低质量、稀疏文库或抽样后 barcode 支持量不足会导致假阴性证据；多组学丰度不同也可能让指标冲突，不能将不一致等同于错配。

## 并行与资源

- `--threads N`：并行扫描的 FASTQ 文件数上限，默认 `min(CPU 核数, 输入文件数, 8)`；`--threads 1` 强制串行。所有样本的 cDNA / oligo 文件一起进入同一个有界线程池，不再按样本排队。
- 并行粒度是**文件**：每个文件独立完成结构检测、全流扫描、蓄水池抽样和局部计数；局部计数即时合并（整数加法，与到达顺序无关），文件统计按原顺序写回报告，最后统一做 `--min-support` 过滤、范数与 Top-N。抽样种子由规范路径决定，所以 `--threads 1` 与 `--threads 8` 的报告逐字节一致；集成测试覆盖全量与抽样两种模式。
- **单个 FASTQ 内部仍串行**：gzip/DEFLATE 是有状态顺序流，蓄水池抽样依赖记录序号，按字节切分或并行抽样会破坏可复现性。想更快可显式 `--sampling prefix`（顺序敏感），或把数据按 lane 拆成多个文件以利用文件级并行。
- 内存随并发数放大：每个在跑的文件各持有一个蓄水池（默认上限下约 8 MB）和一张局部计数表，扫描结果即时合并后释放；同时所有文库的合并表都要常驻以做全交叉比较。文件多、barcode 种类多时建议手动限制 `--threads`，例如 4。
- 实测（Apple M3 / 8 核、8 个 gzip lane 文件共 400 万条 R1、`--max-records-per-file 0`、热缓存）：`--threads 1` 0.81 s / 12 MiB，`--threads 4` 0.29 s / 24 MiB，`--threads 8` 0.23 s / 34 MiB。收益取决于文件数；单个大文件仍受单流 gzip 限制。

## 配对判定

目录结构给出「预期配对」：`<sample>/cDNA` ↔ `<sample>/oligo`。barcode 证据给出「实际领先候选」。二者对照形成每个样本一行：

| `pairings[].assessment` | 含义 |
|---|---|
| `expected_pair_supported` | 三个指标一致领先的同名 oligo，且相关候选均完整计数；仍未校准 |
| `suspected_cross_sample_pairing` | 领先候选属于别的样本，目录约定疑似错位 |
| `unresolved` | 证据不足 / 指标冲突 / 并列 / 抽样未确认，没有唯一领先候选 |
| `expected_candidate_missing` | 该样本没有对应的 oligo 文库 |

`pairings[]` 还记录预期配对的共有 barcode 数、跨样本候选中共有 barcode 最多者（无共有 barcode 时为空；只是计数观察，不是排名结论）。`pairs[]` 中 `expected_pair=true` 表示同一目录内的组合，并带 `cdna_sample` / `oligo_sample`；`rankings[]` 保留 `status`、`preferred_candidate`、`expected_candidate` 与 `expected_pair_assessment`。除 `metrics_agree_unvalidated` 状态外，`preferred_candidate` 一律为空。

## Barcode 提取与独立过滤

| 结构 | CB1（1-based） | CB2（1-based） |
|---|---|---|
| `standard`：无暗反应 | R1 7–16 | R1 23–32 |
| `dark`：暗反应 | R1 1–10 | R1 11–20 |

两边分别提取完整 `CB1 + CB2`，不先使用候选 oligo 校正 cDNA。默认要求 barcode 的每个碱基 Phred+33 ≥20、两段各自精确命中内置白名单；含 N、低质量、过短和白名单外的记录分别统计。不进行错配纠正，避免人为增加交集。

白名单来自本地 `dnbc4tools_py/config/cellbarcode/cell_barcode_1536.txt.gz`，以文本形式随项目保存在 `assets/barcode_whitelist.txt`（1536 段 10 bp barcode）。仅针对这套 C4 结构；**不支持 scRNAv1HT 无暗反应结构、C6、自定义结构或其他白名单**。

默认 `auto` 对每个文件前最多 1 万条记录同时测试两种坐标，采用独立质量过滤及白名单命中数进行判断。最高命中率需 ≥10%，命中数需 ≥`min(20, 探测记录数)`，且至少是另一种结构的 3 倍。证据不足时明确报错，不静默选择结构。这些是结构检测规则，不是样本配对阈值。确认结构后可手工指定：

```bash
--cdna-layout standard --oligo-layout dark
```

## 指标与解释

默认每个文库内完整 barcode 至少有 3 条有效 reads 才参与比较（`--min-support`）；同一样本的多个 lane 支持量会相加。过滤前后的数量和每个输入文件的布局、质量过滤计数都写入 JSON。`scanned_records` / `scanned_valid_records` 对应实际扫描流；`sampled_records` 对应保留的抽样，`valid_records` 与拒绝原因计数对应这份抽样；`input_fully_scanned` 表示是否读到 EOF。支持量过滤作用于抽样后计数，并非全文件分子数。

| JSON 指标 | 定义 |
|---|---|
| `shared_barcodes` | 两边支持量达标的完整 barcode 交集数量 |
| `barcode_jaccard` | 交集数量 / 并集数量 |
| `cdna_barcode_coverage` | 交集数量 / cDNA 达标 barcode 数量 |
| `oligo_barcode_coverage` | 交集数量 / oligo 达标 barcode 数量 |
| `cdna_read_coverage` | 交集 barcode 的 cDNA reads / 全部达标 barcode 的 cDNA reads |
| `oligo_read_coverage` | 相反方向的 reads 覆盖率 |
| `count_cosine` | 两边完整 barcode read 计数向量的余弦相似度 |
| `top_barcode_overlap` | 两边 top barcode 交集 / 两边 top 集合大小的较小值 |

Top 默认各取最多 10000 个 barcode（`--top-barcodes`），相同支持量按 barcode 编码排序。分母为零的指标为 JSON `null`、终端 `NA`。

### 多指标一致性，而非单向覆盖率排名

只对两边各至少有 100 个达标 barcode 的组合比较候选（`--min-distinct-barcodes`）。这个数字仅是数量门槛，不是统计学置信阈值。

v0.2 起分别计算 `cdna_read_coverage`、`barcode_jaccard`、`count_cosine` 的候选顺序和领先候选集合。不再用任意加权分数合并指标。`metric_orders` 是三个各自的顺序，`metric_leaders` 明确保留指标内并列（差值 ≤1e-12 视为数值并列）；`consensus_leaders` 是三组领先集合的交集。

- `metrics_agree_sampled_unconfirmed`：三个指标有唯一共同领先候选，但涉及部分记录抽样或未完整扫描，暂不填充 `preferred_candidate`。
- `metrics_agree_unvalidated`：至少两个候选满足数量门槛、相关候选均完整计数，且三个指标的领先集合只有一个共同候选；填充 `preferred_candidate`，仍不代表已确认配套。
- `conflicting_metrics`、`tied_candidates`、`single_candidate_unresolved`、`insufficient_data`、`no_shared_barcodes`：详见 `rankings[].status`。

共享白名单会在不同实验中重复使用，深度很大的 oligo 可能覆盖大量 barcode，PCR 重复会影响 reads 权重；即使三个指标一致，也可能存在系统性偏差。因此优先候选仍需用真实已知正确/故意错配数据校准。程序不输出 p 值、匹配概率或生物学身份结论。

## 输出与验证

终端输出配对指标 TSV（含 `expected_pair` 列），stderr 输出每个样本的配对结论摘要，JSON 包含 `schema_version: 3`、参数、`samples[]`（发现的目录与 R1/R2 文件）、`cdna[]` / `oligo[]` 文库统计、`pairs[]` 配对矩阵、`rankings[]`、`pairings[]` 和解释限制。完成全部分析后通过临时文件发布 JSON，默认拒绝覆盖已有文件；显式 `--force` 可替换。解析失败时不替换已有报告。

```bash
cargo test --locked
cargo clippy --all-targets --locked -- -D warnings
cargo fmt --check
```

测试覆盖正确/错配候选、跨样本交换、目录别名与多 lane、R2 忽略、低质量、零分母、重复输入、损坏 FASTQ/gzip 和输出保留。输入错误以非零退出；配对不确定是正常报告结果，不会以错误退出。

## 模拟数据演示

以下目录必须尚不存在：

```bash
cargo run --example make_demo -- /tmp/library-pair-demo
cargo run --release -- --input /tmp/library-pair-demo --output /tmp/library-pair-demo/report.json
```

演示生成 `SampleA`、`SampleB`、`SampleC` 三个样本目录，各 120 个完整 barcode、每个 3 条 reads；A/B 的 oligo 文件互相交换，C 正确配对。预期 A、B 为 `suspected_cross_sample_pairing`，C 为 `expected_pair_supported`。这只是实现验证，不代表真实样本的指标分布。

## 版本变化

- v0.3：CLI 从 `--cdna NAME=PATH` / `--oligo NAME=PATH` 改为样本目录 `--input DIR`；JSON `schema_version` 升为 3，新增 `samples[]`、`pairings[]`、`pairs[].expected_pair`、`pairs[].cdna_sample`、`pairs[].oligo_sample`、`rankings[].expected_candidate`、`rankings[].expected_pair_assessment`，以及 `cdna[]` / `oligo[]` 的 `sample` / `role` 字段。旧参数与 schema 2 报告不再兼容。另支持 FASTQ 符号链接，并新增文件级并行的 `--threads`（报告与串行逐字节一致）。
- v0.2：从单向覆盖率排名改为多指标一致性，默认抽样改为全流蓄水池。速度、峰值内存实测见 [`docs/assessment-2026-11-06.md`](docs/assessment-2026-11-06.md)，行为与迁移说明见 [`docs/changes-v0.2.md`](docs/changes-v0.2.md)；这两份文档描述旧版输入模型。

保留旧报告作为结果记录，不应将跨版本的排序差异直接解释为样本变更。全流、不抽样时计数和配对指标应保持一致；余弦范数与点积使用整数累计后转浮点，最后几位可能不同。

解析缓冲区循环复用；Top-N 使用有界堆（O(K) 额外空间）；每文库范数仅计算一次；过滤后在表明显稀疏时收缩 HashMap。蓄水池额外存储每条保留记录的 8 字节 barcode/拒绝代码，默认上限下每文件约 8–8.4 MB 的缓冲区容量，处理完当前文件后释放；不是将完整 FASTQ 读入内存。

独立抽样会丢失稀有 barcode 的支持记录：在每 barcode 仅 3 条 reads 的稀疏数据中，抽样后再按 3 条支持量筛选，可能大幅降低同一完整分布的重合率。部分抽样的指标仅作探索，追加扫描/完整计数后再比较；不要将低覆盖率直接视为错配。

## 基准

`python3 benchmarks/run_macos.py [--output-dir DIR]` 会生成样本目录结构的合成 R1 数据，依赖已构建的 release 二进制与 macOS `/usr/bin/time -l`，结果保留在目录内。脚本固定 `--threads 1`，测的是单线程单文件吞吐，便于与历史数据比较；文件级并行的收益见「并行与资源」。基准脚本测量的是实现性能，不作生物学验证。
