# v0.2 调整与验证（2026-11-06）

## 行为变化

- 不再用单向覆盖率直接决定第一候选。独立展示 cDNA reads 覆盖率、Jaccard、cosine 的排名及领先集合；没有共同领先者时返回 `conflicting_metrics`，不填充优先候选。
- 指标并列、证据数量不足、无交集、单候选、部分记录抽样分别保留为不同状态，不强行选择。即使完整计数的三个指标一致，也只标记未校准的优先候选，不证明技术配套或生物学身份。
- 默认使用固定种子蓄水池抽样：每文件读到 EOF，随机保留最多指定数量记录；不同文库/文件使用独立随机流。同配置可复现，但改变记录顺序会改变具体抽样。`--sampling prefix` 显式选择旧的快速前缀读取；`--max-records-per-file 0` 保留全部记录。
- JSON schema 从 1 升到 2，删除单指标 `candidates` / `top_coverage_ties`。使用 `eligible_candidates`、`metric_orders`、`metric_leaders`、`consensus_leaders`、`preferred_candidate` 并结合 `status`。合格候选名称列表不是优先级顺序。
- JSON 区分扫描量与抽样量；每文件记录实际采样种子及是否读到 EOF。质量过滤计数是对保留记录统计，`scanned_valid_records` 是全扫描流的有效记录数。

## 性能改动

- FASTQ 解析缓冲区复用，直接借用序列/质量字节，不为每条记录复制字符串。
- 有界堆维护 Top-N，O(B log K) 时间、O(K) 额外空间；保留原先相同支持量按 barcode 编码排序的规则。
- 每文库范数只计算一次，范数平方及点积用 u128 累计后转浮点，避免重复扫描与遍历顺序造成的浮点末位变化。
- 支持量过滤后释放空 HashMap，明显稀疏时收缩容量；输入计数阶段的单表峰值仍不可消除。

## 相同数据对比

Apple M3 / macOS，release，每个版本交错重复 3 次；耗时取中位数、RSS 取最大值。全部采用完整计数 `--max-records-per-file 0`，不是默认的限量蓄水池模式。输入为预生成且缓存已热的合成 R1；不能推广到真实质量分布、Linux、冷缓存或网络盘。4×4 使用同内容输入的压力负载，不用于判断真实配对。

| 数据 | 总 R1 记录数 | 耗时秒：v0.1 → v0.2 | 峰值 RSS MiB：v0.1 → v0.2 |
|---|---:|---:|---:|
| small_plain | 200,000 | 0.079 → 0.041 | 4.61 → 4.70 |
| medium_plain | 2,000,000 | 0.651 → 0.310 | 18.77 → 15.80 |
| medium_gzip | 2,000,000 | 0.692 → 0.389 | 18.88 → 15.95 |
| large_gzip | 3,600,000 | 1.385 → 0.742 | 63.33 → 54.98 |
| singletons_gzip | 2,000,000 | 0.932 → 0.524 | 104.98 → 70.81 |
| four_by_four_gzip | 8,000,000 | 2.774 → 1.446 | 45.62 → 42.91 |

原始重复测量见 `benchmarks/v0.2-comparison.json`。逐组合验证完整数据的交集、Jaccard、双向 reads 覆盖率、cosine 和 Top-N 指标差异 <1e-12；候选决策按新规则有意改变。

## 算法边界回归

1. 宽泛候选：cDNA 200 个 barcode，预期候选覆盖 180 个，背景候选覆盖全部 200 个并另含 1800 个。原程序选背景；新程序返回 `conflicting_metrics`，`preferred_candidate=null`。
2. 顺序反例：两文件完整分布相同、前部 barcode 不重合。各取前 300 条、min_support=1 时旧前缀覆盖率 0；当前规范路径与固定种子下，独立蓄水池覆盖率 0.860。重复运行结果一致。随机抽样不保证覆盖率等于完整数据，也不保证重新排列记录后的样本相同。
3. 完全无交集、指标并列、低质量及零分母都不会产生优先候选。
4. 抽样上限之后的损坏尾部：默认蓄水池仍会发现并失败，显式前缀模式只验证读取部分。

### 稀疏数据不能依赖抽样后的绝对覆盖率

在两份各 180 万记录、60 万 barcode（每个恰好 3 条）的相同分布中，各独立抽取 100 万条，再按 min_support=3 过滤，会分别只保留约 10.27 万个 barcode；cDNA 覆盖率约 0.17。完整计数时覆盖率为 1。这是浅抽样与支持量门槛共同造成的真实限制，不是错配证据。

因此部分记录抽样时，即使多个候选的指标一致也返回 `metrics_agree_sampled_unconfirmed`，只展示共同领先者，不填充 `preferred_candidate`。需追加样本量或完整统计后再比较。这个案例的默认蓄水池模式实测中位耗时约 0.666 秒，峰值 RSS 50.30 MiB；它仍扫描全部 360 万条记录，只是减少保留和统计量，不能承诺比读取前缀更快。

## 使用建议

```bash
# 默认：全流读取、每文件蓄水池保留最多 100 万记录
library-pair-check --cdna RNA=cDNA_R1.fq.gz --oligo O=oligo_R1.fq.gz --output preview.json

# 确认：保留全部记录，仍只给未校准的指标证据
library-pair-check --cdna RNA=cDNA_R1.fq.gz --oligo O=oligo_R1.fq.gz --max-records-per-file 0 --output full.json

# 仅作快速预览：显式前缀，可能受记录顺序影响
library-pair-check --cdna RNA=cDNA_R1.fq.gz --oligo O=oligo_R1.fq.gz --sampling prefix --output prefix.json
```
