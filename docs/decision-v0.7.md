# v0.7 丰度选择与集合诊断分离

候选判断只使用 `count_cosine` 和 `log_count_cosine` 的共同领先集合。Jaccard、cDNA reads覆盖率仍保留全部数值、排序及领先集合，但不再否决丰度候选。没有新增扫描、浮点指标或饱和度切换阈值；配对数值计算与v0.6保持一致。

## 状态与输出

- `abundance_agree_unvalidated`：至少两个合格候选、丰度指标有唯一共同领先者、该候选至少100个共有达标barcode，且合格候选均完整扫描，才填充 `preferred_candidate`。
- `abundance_agree_sampled_unconfirmed`：存在未完整扫描的合格候选，保留描述性领先者但不填充优先候选。
- `conflicting_abundance_metrics`：两种cosine的领先集合无交集。
- 共有数量不足、无交集、数量不足、并列与单候选保护继续保留。

新增 `rankings[].abundance_candidate` 表示描述性唯一丰度领先者，弱交集和抽样情况下也可能非空，不能代替 `preferred_candidate`。只有一个合格候选时为空。

`overlap_candidate` 是Jaccard的唯一领先者，并列则为空；`overlap_agrees` / `coverage_agrees` 检查丰度领先者是否包含在相应诊断指标的领先集合中。并列包含时为true，表示兼容；没有唯一丰度领先者时为null，而非false。诊断警告为 `jaccard_disagrees_with_abundance` / `read_coverage_disagrees_with_abundance`，同步到 `pairings[]` 与终端摘要。

`consensus_leaders` 现在仅表示两种cosine的交集，不再是四指标交集。`informative_metrics` 仍保留四指标区分力说明；`ranking_basis` 只描述丰度选择依据，为 `abundance_only`、`no_discriminating_abundance_metrics` 或 `descriptive_only`。

两种cosine相关，不是独立验证。没有已校准的绝对分值或最小领先幅度阈值；`expected_pair_supported` 仅表示未校准的丰度证据支持目录配对。单对仍为 `single_candidate_unresolved`，不输出匹配概率。

## 已提供真实报告的复核

使用提供的v0.6报告中121组已保存pair指标建立回归fixture，样本标签替换为deep_11至deep_21，不包含原始FASTQ路径。用生产代码实际执行新排名函数，逐项检查：

- 11个cDNA均选择同号oligo，状态 `abundance_agree_unvalidated`。
- 14、16、17、18、20、21的Jaccard领先者仍为13号，`overlap_agrees=false`，保留明确警告。
- 其余5组 `overlap_agrees=true`；全部11组 `coverage_agrees=true`。

这是对保存指标的判定回放，没有重新读取55.39亿条FASTQ，也没有把目录配对当作生物学真值来校准阈值。

## 验证与兼容

19项单元测试和12项CLI集成测试通过，clippy、fmt与release构建通过。除真实指标回归外，覆盖宽泛背景不会因覆盖率高而自动胜出、诊断指标无法打破丰度并列、诊断并列的兼容语义、终端警告、两种cosine冲突、弱交集、部分扫描、单对以及串并行一致。

工具版本0.7.0，JSON schema升级为7。消费者需更新状态名、`consensus_leaders`与`ranking_basis`解释；终端摘要新增诊断标记，TSV数值列不变。历史v0.6文档仅描述旧规则，应以当前README为准。
