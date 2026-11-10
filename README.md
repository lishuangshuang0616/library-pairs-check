# library-pair-check

用于 **DNBelab C4 多样本 cDNA / oligo 原始 R1 FASTQ** 的 barcode 指纹配对检查，Rust v0.7。

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

## 运行

```bash
./target/release/library-pair-check \
  --input /data/project \
  --output pairing-report.json
```

默认**统计每个 R1 文件的全部 reads**。若只想快速预览，用 `--reads N` 让每个文件只读前 N 条：

```bash
./target/release/library-pair-check \
  --input /data/project --reads 1000000 \
  --output preview.json
```

`--reads N` 会明显加快大文件的分析，但只读前缀：对记录顺序敏感，未读尾部不做 gzip 完整性 / FASTQ 解析检查；因为没读到 EOF，数量充足且两种 cosine 一致时排名状态会是 `abundance_agree_sampled_unconfirmed`，**不填充 `preferred_candidate`**（领先者仍可在 `rankings[].consensus_leaders` 查看）。`--reads 0` 与不传该参数等价。N 也不宜太小：内置要求每个 barcode 至少 3 条 reads、每库至少 100 个达标 barcode，N 只有几百条时会出现 `insufficient_data`，一般建议 N ≥ 1000。

只有一对输入时，数量充足且有交集也会返回 `single_candidate_unresolved`；程序不会仅因唯一候选而确认配套。

所有 cDNA × 所有 oligo 都会比较，不做强制一对一分配。每个文件独立检测 barcode 结构、独立执行 reads 上限；同一文件内部应保持结构一致。低质量、稀疏文库或 barcode 支持量不足会导致假阴性证据；多组学丰度不同也可能让指标冲突，不能将不一致等同于错配。

## 参数

| 参数 | 默认 | 说明 |
|---|---|---|
| `--input DIR` | 必填，可重复 | 样本根目录或单个样本目录 |
| `--output PATH` | 必填 | JSON 报告；**已存在则直接覆盖**（原子写，分析失败不会破坏旧报告） |
| `--reads N` | `0` = 全部 | 每个 R1 FASTQ 只读前 N 条；大于 0 时顺序敏感、尾部不校验 |
| `--threads N` | `min(CPU 核数, 文件数, 8)` | 并行扫描的 FASTQ 文件数上限；`1` 为串行 |
| `--min-support N` | `3` | 完整 barcode 至少多少条 reads 才参与比较 |
| `--cdna-layout` / `--oligo-layout` | `auto` | `auto` / `standard` / `dark`，见「Barcode 提取」 |

以下阈值**内置固定、不可配置**：barcode 每个碱基最低 Phred+33 = 20（Q20）、Top barcode 集合 10000 个、参与排名的最低达标 barcode 数 100、1536 段 10 bp 白名单。

## 并行与资源

- 并行粒度是**文件**：所有样本的 cDNA / oligo 文件一起进入同一个有界线程池，不再按样本排队。每个文件独立完成结构检测、扫描和局部计数，局部计数通过零容量有界通道交给主线程合并（每个工作线程最多保留一份局部表，首份表直接接管；整数加法，与到达顺序无关），文件统计按原顺序写回报告，最后统一做 `--min-support` 过滤、范数与 Top-N。因此 `--threads 1` 与多线程的报告（除 `timings` 外）完全一致，集成测试覆盖全量与 `--reads N` 两种模式。
- **单个 FASTQ 内部仍串行**：gzip/DEFLATE 是有状态顺序流，不能按字节切给多线程。想更快可用 `--reads N`，或把超大文件按记录块拆成多个文件以利用文件级并行。
- 内存来自每个在跑文件的局部计数表和所有文库的合并表：随并发数、样本数、观察到的不同 barcode 数增长。barcode 空间上限是 \(1536^2 \approx 236\) 万种，超大文件全量计数时计数表可能明显大于抽样时代；内存吃紧就降低 `--threads`，或用 `--reads N` 限制输入量。报告中给出的 `timings` 可帮助定位耗时来源。
- 实测（Apple M3 / 8 核、8 个 gzip lane 文件共 400 万条 R1、全量、热缓存）：`--threads 1` 0.81 s / 12 MiB，`--threads 4` 0.29 s / 24 MiB，`--threads 8` 0.23 s / 34 MiB（v0.3 抽样实现的数字，量级仍可参考）。收益取决于文件数；单个大文件仍受单流 gzip 限制。

## 阶段计时

每次运行会在 stderr 打印一行阶段耗时摘要和最慢的 5 个文件，JSON 顶层也写入 `timings`：

| 字段 | 含义 |
|---|---|
| `discover_seconds` | 样本目录与 FASTQ 发现 |
| `scan_wall_seconds` | 全部扫描、局部合并与 profile 组装的墙钟时间（包含 `merge_seconds`，不可相加） |
| `merge_seconds` | 局部计数合并、支持量过滤、范数与 Top-N；是 `scan_wall_seconds` 的重叠子集 |
| `compare_seconds` | 所有 cDNA × oligo 的全交叉指标计算 |
| `rank_seconds` | 每个 cDNA 的多指标排名与 `pairings` 组装 |
| `total_seconds` | 从启动到构造报告前，不含 JSON 序列化与落盘 |
| `threads` / `scanned_records` / `records_per_second` | 线程数、总扫描 R1 记录数、整体吞吐 |

`timings.files[]` 是每个文件的细分：`detect_seconds`（结构检测）、`read_seconds`（打开、解压、解析、barcode 提取与计数）、`scan_seconds`（该文件总耗时）、`scanned_records`。并行时各文件 `scan_seconds` 会重叠，其和大于 `scan_wall_seconds`。

判断瓶颈：

- `scan_wall_seconds` 占绝对大头：用 `pigz -dc <R1> | wc -l` 或 `cat <R1> > /dev/null` 测裸读。裸读就慢说明瓶颈在共享存储 I/O 或单流 gzip 解压，加线程无用；裸读明显更快则是解析开销。
- `compare_seconds` 大：样本数或每库 barcode 种类多，瓶颈在全交叉比较（目前串行）。
- `records_per_second` 与本机单线程基线（合成数据约 5–6 M records/s）对比，可判断是 CPU 还是 I/O 限制；文件很小时固定开销会拉低该值，不必解读。

`timings` 只是诊断数据，不参与结果一致性：除该字段外，`--threads 1` 与多线程的报告完全一致。

## 配对判定

目录结构给出「预期配对」：`<sample>/cDNA` ↔ `<sample>/oligo`。barcode 证据给出「实际领先候选」。二者对照形成每个样本一行：

| `pairings[].assessment` | 含义 |
|---|---|
| `expected_pair_supported` | 两种 cosine 有唯一共同领先的同名 oligo，且相关候选均完整扫描；仍未校准 |
| `suspected_cross_sample_pairing` | 领先候选属于别的样本，目录约定疑似错位 |
| `unresolved` | 证据不足 / 指标冲突 / 并列 / 未完整扫描，没有唯一领先候选 |
| `expected_candidate_missing` | 该样本没有对应的 oligo 文库 |

`pairings[]` 还记录预期配对的共有 barcode 数、跨样本候选中共有 barcode 最多者（无共有 barcode 时为空；只是计数观察，不是排名结论）。`pairs[]` 中 `expected_pair=true` 表示同一目录内的组合，并带 `cdna_sample` / `oligo_sample`；`rankings[]` 保留 `status`、`preferred_candidate`、`expected_candidate` 与 `expected_pair_assessment`。除 `abundance_agree_unvalidated` 状态外，`preferred_candidate` 一律为空。

## Barcode 提取与过滤

| 结构 | CB1（1-based） | CB2（1-based） |
|---|---|---|
| `standard`：无暗反应 | R1 7–16 | R1 23–32 |
| `dark`：暗反应 | R1 1–10 | R1 11–20 |

两边分别提取完整 `CB1 + CB2`，不先使用候选 oligo 校正 cDNA。固定要求 barcode 的每个碱基 Phred+33 ≥20、两段各自精确命中内置白名单；含 N、低质量、过短和白名单外的记录分别统计。不进行错配纠正，避免人为增加交集。

白名单来自本地 `dnbc4tools_py/config/cellbarcode/cell_barcode_1536.txt.gz`，以文本形式随项目保存在 `assets/barcode_whitelist.txt`（1536 段 10 bp barcode）。仅针对这套 C4 结构；**不支持 scRNAv1HT 无暗反应结构、C6、自定义结构或其他白名单**。

默认 `auto` 对每个文件前最多 1 万条记录同时测试两种坐标。最高命中率需 ≥10%，命中数需 ≥`min(20, 探测记录数)`，且至少是另一种结构的 3 倍。证据不足时明确报错，不静默选择结构。这些是结构检测规则，不是样本配对阈值。确认结构后可手工指定：

```bash
--cdna-layout standard --oligo-layout dark
```

## 指标与解释

默认每个文库内完整 barcode 至少有 3 条有效 reads 才参与比较（`--min-support`）；同一样本的多个 lane 支持量会相加。每个输入文件的布局、扫描记录数、是否读到 EOF、有效记录数与各类拒绝计数都写入 JSON。`input_fully_scanned=false` 表示被 `--reads N` 截断，此时结论视为未确认。

| JSON 指标 | 定义 |
|---|---|
| `shared_barcodes` | 两边支持量达标的完整 barcode 交集数量 |
| `barcode_jaccard` | 交集数量 / 并集数量 |
| `cdna_barcode_coverage` | 交集数量 / cDNA 达标 barcode 数量 |
| `oligo_barcode_coverage` | 交集数量 / oligo 达标 barcode 数量 |
| `cdna_read_coverage` | 交集 barcode 的 cDNA reads / 全部达标 barcode 的 cDNA reads |
| `oligo_read_coverage` | 相反方向的 reads 覆盖率 |
| `count_cosine` | 两边完整 barcode read 计数向量的余弦相似度；整体倍数缩放不改变该值 |
| `log_count_cosine` | 对达标 barcode 计数作 `ln(1+count)` 后的余弦；缺失 barcode 为 0，降低极高计数支配，但不具备深度缩放不变性 |
| `abundance_overlaps[]` | Top 100 / 1000 / 10000 的交集除以较小集合大小；边界同分全部纳入，附实际集合大小 |
| `top_barcode_overlap` | 两边 top barcode 交集 / 两边 top 集合大小的较小值 |

Top 固定各取最多 10000 个 barcode，相同支持量按 barcode 编码排序。分母为零的指标为 JSON `null`、终端 `NA`。

### 丰度决定候选，集合指标单独诊断

只对两边各至少有 100 个达标 barcode 的组合排名；所有四个指标都保留各自的 `metric_orders` 和 `metric_leaders`，但只有 `count_cosine` 与 `log_count_cosine` 参与候选选择。`consensus_leaders` 是这两组领先集合的交集，差值 ≤1e-12 视为并列。两种 cosine 是相关证据，不能当作两项独立验证。

- `abundance_agree_unvalidated`：至少两个合格候选，两种 cosine 有唯一共同领先者，其共有达标 barcode 至少100，且所有合格候选完整扫描；填充 `preferred_candidate`，仍未校准。
- `abundance_agree_sampled_unconfirmed`：同样有唯一丰度领先者、共有数量达标，但存在未完整扫描的合格候选；不填充优先候选。
- `conflicting_abundance_metrics`：两种 cosine 没有共同领先候选，不选择。
- `insufficient_shared_barcodes`：唯一丰度领先者的共有 barcode 不足100；数量保护，不是匹配阈值。
- `tied_candidates`、`single_candidate_unresolved`、`insufficient_data`、`no_shared_barcodes` 保留原有边界。

`abundance_candidate` 是至少两个合格候选时的唯一共同丰度领先者；弱交集或抽样时仍可作为描述性线索，不能替代 `preferred_candidate`。单候选时该字段为空。`overlap_candidate` 是 Jaccard 的唯一领先者，并列时为空。

`overlap_agrees` / `coverage_agrees` 表示丰度领先者是否在 Jaccard / cDNA reads覆盖率的领先集合中。诊断指标并列且包含它时为 true，表示兼容，不代表唯一支持；没有唯一丰度领先者时为 null。诊断冲突写入 `diagnostic_warnings`：`jaccard_disagrees_with_abundance` / `read_coverage_disagrees_with_abundance`。这些字段同时出现在 `pairings[]` 和终端摘要中，Jaccard与覆盖率不否决丰度候选。

`informative_metrics` 继续描述四指标中哪些能区分候选，只有两种cosine决定选择；`ranking_basis` 为 `abundance_only`、`no_discriminating_abundance_metrics` 或 `descriptive_only`（少于两个合格候选）。`count_cosine_margin` / `log_count_cosine_margin` 是各自前两名差值；没有校准的最小差值阈值。`coverage_margin` 仅作诊断。

深测序时cDNA可能覆盖绝大部分barcode空间，Jaccard会受oligo集合大小影响。程序不按某个占用率阈值切换规则，始终分别报告丰度判断与集合诊断。barcode空间占用率不等于常规测序饱和度。

log cosine不具备整体深度倍数缩放不变性；大量barcode存在且丰度平坦时，两种cosine都为1也不能确认配套。`count_effective_barcodes` 为 `(sum counts)^2 / sum(count²)`，越接近达标barcode总数表示越平坦；`largest_barcode_read_fraction` 显示最大barcode的reads占比。Top分层只作诊断，不参与投票。

只有一对输入时仍保持 `single_candidate_unresolved`、`preferred_candidate=null`。程序不输出p值、匹配概率或生物学身份结论，仍需真实正确与故意错配数据校准；共享白名单、PCR和深度都可能造成系统偏差。

## 输出与验证

终端输出配对指标 TSV（含 `expected_pair` 列），stderr 输出每个样本的配对结论摘要、阶段耗时和最慢文件，JSON 包含 `schema_version: 7`、`settings`、`samples[]`（发现的目录与 R1/R2 文件）、`cdna[]` / `oligo[]` 文库统计、`pairs[]` 配对矩阵、`rankings[]`、`pairings[]`、`timings` 和解释限制。完成全部分析后通过临时文件原子发布 JSON，**默认覆盖同名输出**；解析失败时不替换已有报告。

```bash
cargo test --locked
cargo clippy --all-targets --locked -- -D warnings
cargo fmt --check
```

测试覆盖正确/错配候选、跨样本交换、`--reads` 截断、目录别名与多 lane、R2 忽略、低质量、零分母、重复输入、损坏 FASTQ/gzip、旧参数拒绝与输出覆盖。输入错误以非零退出；配对不确定是正常报告结果，不会以错误退出。

## 模拟数据演示

以下目录必须尚不存在：

```bash
cargo run --example make_demo -- /tmp/library-pair-demo
cargo run --release -- --input /tmp/library-pair-demo --output /tmp/library-pair-demo/report.json
```

演示生成 `SampleA`、`SampleB`、`SampleC` 三个样本目录，各 120 个完整 barcode、每个 3 条 reads；A/B 的 oligo 文件互相交换，C 正确配对。预期 A、B 为 `suspected_cross_sample_pairing`，C 为 `expected_pair_supported`。这只是实现验证，不代表真实样本的指标分布。

## 版本变化

- v0.7：候选选择仅由原始与log cosine共同决定，Jaccard和reads覆盖率改为诊断。新增丰度候选、集合候选、一致性字段和终端告警；状态改为 `abundance_agree_unvalidated` / `abundance_agree_sampled_unconfirmed` / `conflicting_abundance_metrics`。JSON schema升为7，所有配对数值计算保持不变。验证见 [`docs/decision-v0.7.md`](docs/decision-v0.7.md)。

- v0.6：保留原始 cosine，新增 log-count cosine、包含边界并列的 Top 100/1000/10000 重合、丰度集中度与候选区分依据；log cosine 加入共同领先规则。JSON schema 升为 6，TSV 新增 `log_count_cosine` 列。单对仍只报告证据；指标冲突或数据不足不确认。验证见 [`docs/abundance-v0.6.md`](docs/abundance-v0.6.md)。

- v0.5：有界结果通道与首表接管降低内存峰值；配对遍历较小计数表，排名与摘要只扫描对应 cDNA 的候选行；增加共有 barcode 至少 100 的数量保护和 `insufficient_shared_barcodes` 状态。JSON schema 升为 5，新增 `pairs[].sufficient_shared_barcode_count`。单对仍为 `single_candidate_unresolved`，没有新增匹配概率或已校准阈值。优化实测见 [`docs/optimization-v0.5.md`](docs/optimization-v0.5.md)。

- v0.4：CLI 大幅简化。删除蓄水池抽样与 `--sampling` / `--max-records-per-file` / `--seed`，改为单一 `--reads N`（默认 `0` = 全部 reads，否则每个文件只读前 N 条）；删除 `--min-base-quality`（固定 Q20）、`--top-barcodes`（固定 10000）、`--min-distinct-barcodes`（固定 100）、`--force`（默认覆盖输出）。JSON `schema_version` 升为 4，`settings` 只剩 `reads_per_file` 与 `min_support`，`files[]` 去掉 `sampled_records` / `sampling_seed` / `scanned_valid_records`，`timings.files[]` 去掉 `accumulate_seconds`。旧参数与 schema 3 报告不兼容。
- v0.3：CLI 从 `--cdna NAME=PATH` / `--oligo NAME=PATH` 改为样本目录 `--input DIR`；JSON `schema_version` 升为 3，新增 `samples[]`、`pairings[]`、`pairs[].expected_pair`、`rankings[].expected_candidate` 等字段。另支持 FASTQ 符号链接、文件级并行的 `--threads` 与 `timings` 阶段计时。
- v0.2：从单向覆盖率排名改为多指标一致性，默认抽样改为全流蓄水池。速度、峰值内存实测见 [`docs/assessment-2026-11-06.md`](docs/assessment-2026-11-06.md)，行为与迁移说明见 [`docs/changes-v0.2.md`](docs/changes-v0.2.md)；这两份文档描述旧版输入模型。

保留旧报告作为结果记录，不应将跨版本的排序差异直接解释为样本变更。余弦范数与点积使用整数累计后转浮点，最后几位可能不同。

解析缓冲区循环复用；Top-N 使用有界堆（O(K) 额外空间）；每文库范数仅计算一次；过滤后在表明显稀疏时收缩 HashMap。不是将完整 FASTQ 读入内存。

## 基准

`python3 benchmarks/run_macos.py [--output-dir DIR]` 会生成样本目录结构的合成 R1 数据，依赖已构建的 release 二进制与 macOS `/usr/bin/time -l`，结果保留在目录内。脚本固定 `--threads 1`、使用默认全量 reads，测的是单线程单文件吞吐，便于与历史数据比较；文件级并行的收益见「并行与资源」。基准脚本测量的是实现性能，不作生物学验证。
