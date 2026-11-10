# v0.5 优化与验证

在 2026-11-10 当前 v0.4 工作区基础上修改；保留已有未提交更改。

## 判定行为

保留三个指标的独立排序与共同领先集合。唯一共同领先者的共有达标 barcode 少于 100 时，新增 `insufficient_shared_barcodes`，`preferred_candidate=null`，配对评估为 `unresolved`。`pairs[].sufficient_shared_barcode_count` 展示该数量条件。未将弱候选从排序中删除，避免通过删除竞争者人为制造一致领先。

100 沿用原有每库最低 barcode 数量要求，是最低证据量保护，不是生物学或统计学匹配阈值。即使共享 100 个也可能只是背景。仍需真实正负对照校准；未加入 p 值、匹配概率或单对自动确认。单候选有交集时仍为 `single_candidate_unresolved`。

JSON schema 升为 5，工具版本为 0.5.0；消费者须处理新增字段和状态。原来的弱交集领先结果可能改为无法判定。

## 实现

- 并行结果使用零容量同步通道：不保留排队结果；阻塞发送的线程仍持有自己的局部计数表，因此局部表数量随工作线程数受限。所有文库合并表、正在合并的表及分配器开销仍占内存，线程数不是 RSS 的固定上限。
- 文库首个局部计数表直接接管，避免重新分配和重新插入全表；后续 lane 仍使用整数累加，支持量在合并后过滤。
- 交集遍历较小计数表，并保持 cDNA/oligo reads 权重方向不变。
- 利用配对矩阵的行顺序，每个 cDNA 的排名和摘要只扫描本行，避免反复扫描整个矩阵；全交叉比较仍存在，未引入近似算法。
- 明确 `scan_wall_seconds` 包含合并和 profile 组装；`merge_seconds` 是重叠子集，不能相加。

## 实测

Apple M3 / 16 GiB，macOS release，合成数据热缓存；优化前后交错重复 3 次，外部进程墙钟取中位数，RSS 取最大值。请求 8 线程的单对只有两个文件，实际只有两线程。多样本使用相同内容硬链接，是实现负载测试，不是配对有效性验证。

| 负载 | 线程 | 秒：前 → 后 | MiB：前 → 后 |
|---|---:|---:|---:|
| 200万R1、每库20万barcode | 1 | 0.392 → 0.367 | 21.02 → 15.66 |
| 360万R1、每库60万barcode | 1 | 0.846 → 0.785 | 70.98 → 54.02 |
| 200万R1、每库100万singleton | 1 | 0.725 → 0.575 | 140.05 → 104.94 |
| 360万R1、每库60万barcode | 2 | 0.569 → 0.477 | 88.27 → 69.80 |
| 单样本8个lane文件、800万R1 | 8 | 0.365 → 0.349 | 70.84 → 67.09 |
| 4样本、16配对、80万R1 | 8 | 0.066 → 0.063 | 11.75 → 9.97 |
| 16样本、256配对、320万R1 | 8 | 0.365 → 0.373 | 29.53 → 28.19 |

单线程高复杂度场景峰值内存减少约24–25%。16样本场景未显示明确的耗时改善，不宣称多样本整体提速。合成数据压缩率高，不能外推真实文件、冷缓存或网络磁盘耗时。

所有基准的交集、覆盖率、Jaccard、cosine、Top overlap 及其余既有 pair 字段逐项精确一致。新增边界测试覆盖共享1/99/100、不同大小集合的方向性读数，以及单对完全重合仍不确认。13项单元测试和10项CLI集成测试通过，clippy与fmt检查通过。

## 复现

先保留优化前 release 二进制，再构建当前版本。运行 `benchmarks/run_macos.py` 生成输入；然后执行：

```bash
python3 benchmarks/compare_optimization.py \
  --before /path/to/previous/library-pair-check \
  --after target/release/library-pair-check \
  --data-root /path/to/generated/benchmark-directory \
  --output benchmarks/optimization-v0.5-results.json
```

原始三轮记录保存在 `benchmarks/optimization-v0.5-results.json`。脚本额外生成真正的4/16样本目录，保留临时输入和报告供复核；不删除用户数据。
