# 设计草案:把 quant 方法选择纳入 agentic 优化

> 状态:**设计草案(本轮不实现代码)**。这是一份本地技术文档,供内部评审后再决定是否推进架构改动。对应 issue #48 Point 3。

## 1. 背景与动机

Mokume Plugin 的 MCP service 目前**严格在“蛋白矩阵之后”**工作:它接收一张已经建好的蛋白强度矩阵(`rust/python/mokume/agentic/runner.py`),只在其上调 归一化 / 插补 / DE 方法 / FDR / 效应量门。**quant 方法**(把肽段/特征汇总成蛋白的算法:maxlfq / directlfq / pibaq / topn / sum …)是**上游 `features2proteins` 用 `--quant-method` 定死的**(默认 maxlfq),当前工具从不触碰。

issue #48 指出:quant 方法很可能是**比任何 DE 方法选择影响都大**的杠杆,值得纳入候选空间。维护者原文把全部 ~10 种都列为候选、并说 quant "in particular seems worth promoting",但明确因其是**架构级改动**而"只提出、不直接 PR"。

## 2. 已有证据

- **maxlfq vs directlfq(已实测)**:15 个 spike-in 数据集的评分网格都对比过这两种正经 LFQ 定量算法。示例:PXD000279 maxlfq sensitivity 0.60 vs directlfq 0.49——**quant 轴确实能差出一大截**,支持"值得纳入"。
- **piBAQ / topn / sum(已实测,结论为否)**:后续用 `features2proteins` 从 `quantms.feature.parquet` 在 5 个 LFQ 数据集上**真建了** piBAQ/topn/sum 蛋白矩阵,并跑了同一张网格(limma/deqms × {none, median} × 8 种插补)。这些历史运行使用旧 CLI 名称 `ibaq`,但执行的是现已规范命名为 piBAQ 的 shared-peptide allocation 实现,并非传统 iBAQ。用**分母无关**的度量对比(AUC 与 emp_FDP ≤ 0.05 下的绝对 TP;sensitivity 的分母 gt∩tested 会随 quant/插补 漂移,不可比):

  | 数据集 | best AUC:directlfq / maxlfq | best AUC:piBAQ / topn / sum | 差距 |
  | --- | --- | --- | --- |
  | PXD028735 | 0.929 / 0.953 | 0.671 / 0.560 / 0.547 | +0.28 |
  | PXD007683 | 0.973 / 0.991 | 0.522 / 0.541 / 0.476 | +0.45 |
  | PXD000279 | 0.809 / 0.920 | 0.700 / 0.756 / 0.732 | +0.16 |
  | PXD070049 | 0.967 / 0.972 | 0.916 / 0.915 / 0.935 | +0.04 |

  受控 TP(emp_FDP ≤ 0.05)更悬殊:piBAQ/topn/sum 在 PXD028735 / PXD007683 / PXD000279 **全部归零**,只有 PXD070049 有发现(1280 / 936 / 1049)但仍不及 directlfq 的 1540。**piBAQ/topn/sum 从不胜出。**

- **TMT 上 piBAQ/topn/sum 实测 N/A**:它们是 LFQ 汇总法,`features2proteins` 的 TMT 路径不走这些方法,在 TMT reporter parquet 上产 **0 蛋白**。TMT 的 quant 扫描只对 directlfq/maxlfq 有意义。

- **结论**:现有 **directlfq/maxlfq 双 quant 空间有数据支撑**,**不建议**扩到 piBAQ/topn/sum。

**证据的诚实边界**:

1. 本轮 quant 网格限于 **limma / deqms**;慢的 rots / limrots 为控制时长排除(已记录)。差距量级(AUC 0.04–0.45、TP 归零)远大于换 DE 方法的预期收益,但严格说这三种未在 rots 家族下验证。
2. **PXD003881 为空数据集**(3x/1x 富集差过弱,4v4 下所有方法 0 受控发现),不参与判定。它在 AUC 表上 piBAQ/topn/sum 看似略高(0.645/0.682/0.655 vs 0.530/0.543)属噪声,不构成反证。

## 3. 为什么是"架构级"改动

quant 方法**生产**蛋白矩阵。要在候选空间里比较不同 quant,循环就不能再从"一张现成矩阵"出发,必须**退到肽段/特征层、每个候选重跑定量**。这动到三处:

1. **`CandidateConfig`(`rust/python/mokume/agentic/state.py`)** 增 `quant_method` 字段(str,默认继承上游)。
2. **`runner.run_experiment`** 从"吃现成 `protein_df`"改为"能从特征输入重跑定量":当候选的 `quant_method` 与基线不同时,调用 `features2proteins` 的定量阶段(`QuantificationStage`)重建矩阵,再进入既有的 归一化→插补→DE。
3. **缓存(`PreprocessCache`)** 改造:现在缓存键是 `(normalization, imputation)`,前提是矩阵固定;纳入 quant 后需把 `quant_method` 也纳入缓存键,且缓存对象从"插补后矩阵"扩展到"定量后矩阵"。定量比 归一化/插补 重得多,缓存收益更大但也更占内存。

此外,MCP tool contract 要从 `protein_df` 扩展到“特征层输入 + SDRF(+ FASTA)”。`service.py`、skill contract 和 Policy 相应修改。

## 4. 建议的最小可行设计

- **候选空间**:`quant_method` 作为一等参数;但**不暴力跑全网格**(issue 明确排除)。按第 2 节的实测证据,候选集就是 **maxlfq / directlfq 两种**——piBAQ/topn/sum 已测且从不胜出,不纳入。当前 knowledge graph 会保留完整 quant provenance，但 protein-matrix 入口的 Generation Scope 明确把 quant 标成 frozen axis；真正搜索 quant 仍需本设计所述的特征层输入改造。
- **重定量**:`run_experiment` 增一条"若候选 quant≠已缓存 quant 则重跑 `QuantificationStage`"的分支;否则命中缓存(与现状等价,零回归)。
- **缓存**:键 `(quant, normalization, imputation)`;定量结果单独缓存,内存上限可配。
- **成本控制**:重定量昂贵,建议默认只在**首轮**跨 quant 探索、后续轮固定胜出 quant 只调下游;或提供 `explore_quant: bool` 开关。

## 5. 风险与开放问题(留给评审)

- **计算成本**:每候选重定量,循环变慢一个量级(尤其 maxlfq)。需成本-收益权衡,可能只对小数据集或首轮开启。这是当前**最主要**的开放问题——候选集只剩两种 quant,收益上限也就相应有限。
- **piBAQ 的可比性(已有答案)**:piBAQ 是绝对丰度、sum/topn 是更粗的汇总,第 2 节的实测印证了担心——它们对**组间比值**的 DE 基准确实不合适。"更多 quant = 更好"已被数据否定。
- **纳入范围已由数据定**:maxlfq / directlfq 两种,不是"一次全纳入"。
- **与上游 quantms 的边界**:定量本身在 quantms/mokume 定量层;agentic 只是"选哪个",不重造定量算法。

## 6. 修饰过滤:定位为稳健性维度,不是又一条搜索轴

`features2proteins` 侧还有一个相邻的旋钮:**是否过滤掉带修饰的肽段**(只留非修饰肽段建蛋白矩阵)。它和 quant 一样发生在"蛋白矩阵之前",容易被顺手当成候选空间的又一个轴。**不应该这么做。**

- **不是搜索轴**:如果把"去修饰 / 不去修饰"塞进候选空间、按分数择优,那就是在**用一个和处理组无关的技术开关去最大化 n_de**。去修饰会换掉一批参与汇总的肽段,从而抖动矩阵、抖动分数;择优挑中的很可能只是这次抖动的运气方向。这正是 p-hacking 的形状,和"经验 FDR 校准"的目标相反。
- **是稳健性检查**:正确的用法是**先定住配置,再看结论变不变**——同一条胜出流水线,分别在"全部肽段"和"仅非修饰肽段"上跑一遍,比较结论(显著蛋白集合、方向、效应量)是否一致。
  - **不变** → 给结论**加信**:它不是某类修饰肽段的汇总假象。
  - **变了** → **不**自动选分数高的那个,而是在**报告层提示**:该数据集的结论对修饰肽段敏感,列出仅在一侧显著的蛋白,交人判断(可能是真实的修饰驱动生物学,也可能是定位错误 / 定量不稳)。
- **评价指标**:稳健性维度报"一致性"(如两侧显著集合的重叠率、效应量相关),**不报** "哪边 n_de 更多"。n_de 更多在这里不是好消息。
- **成本**:它使流水线跑两遍而非搜索空间翻倍——是**固定 2 倍**开销,不与候选数相乘。可作 `robustness_checks` 开关,默认关闭,在最终报告阶段对**胜出配置**跑一次即可,不进优化循环。

## 7. 归一化轴:同一个架构障碍,外加一个"配置≠事实"的坑

issue #48 Point 3 除 quant 外还点名了归一化搜索空间。`CandidateConfig.normalization` 是**单个字段**,`runner._apply_normalization` 只认 6 个值、且全部作用在**已建好的蛋白矩阵**上。mokume 实际有**三层**归一化:

| 层 | 作用位置 | CLI | agentic 能否触及 |
| --- | --- | --- | --- |
| run | 特征层,跨技术重复/馏分 | `features2proteins --run-normalization`(`FeatureNormalizationMethod`,7 个) | **不能** |
| sample | 肽段层,样本之间 | `features2proteins --sample-normalization`(`PeptideNormalizationMethod`,10 个) | **不能** |
| matrix | 蛋白矩阵层 | `CandidateConfig.normalization`(6 个) | 只搜这层 |

**run 层与 quant 是同一个障碍,不是另一个问题。** 建矩阵时技术重复已被合并,run 这一维在 agentic 的输入里不复存在——不是"没搜",是从这个入口**够不到**。要搜它,就得做第 3 节那三处改动 + 输入契约改造,与 quant 完全同一批。因此 run 归一化**并入本设计**,不单开一份。

sample 层与 matrix 层则做同一件事(拉平样本之间)、只是层级不同,所以 agentic 的 matrix 归一化实际是在**重做**上游 sample 归一化。这也解释了当前 knowledge graph 中 OpDEA 与 Grid evidence 为何常把 `none` 放在参考配置里——那些证据来自“归一化只做一次”的流水线，搬到“上游已归一化过”的矩阵上，`none` 往往就是对的。

### 7.1 上游归一化参数现在是可执行契约,但仍不是 provenance 边车

本节最初记录的“配置接受但计算不读取”问题已经在命令契约修复中消除。当前行为是:

| quant 分支 | `--run-normalization` / `--sample-normalization` 行为 |
| --- | --- |
| `directlfq` | 省略时解析为 `none/none`; 显式活动归一化被拒绝,由 DirectLFQ 自己完成内部归一化 |
| `ratio` | 省略时解析为 `none/none`; 显式活动归一化被拒绝,由每 plex reference ratio 完成缩放 |
| `maxlfq` | 无 dataset-level 方法时走 DirectLFQ-aligned 路径; 请求 quantile 时切换到 built-in MaxLFQ 并实际应用 quantile; 当前不支持的 dataset-level 组合直接报错 |
| `pibaq` | quantile 会实际应用; 其他 dataset-level 方法直接报错 |
| `top<N>` / `sum` 等 cell-based 方法 | 请求的已支持 run/sample 方法均进入实际计算 |

CLI 使用“是否显式传入”的参数来源信息区分默认值和用户请求; API 层也在运行前检查方法作用域。因此不再通过 warning 接受无效配置,也不存在 `--threads` 被 `--directlfq-cores` 覆盖的优先级。

此外 run 归一化在无技术重复时是空操作(`stages.py:440`、`:587`、`:683` 三处都以 `technical_repetitions > 1` 为条件),即对许多数据集,"被冻结的 run-norm 轴"根本不是一条真轴。

**结论:当前 config 已是经过作用域验证的执行意图,矩阵仍是最终计算事实。** 如果未来把 config 序列化为 provenance,仍应记录 resolved/effective 值与实际路由;本轮按既定决定不新增 provenance 边车。

本轮落地的是不依赖上游猜测、对**任何**来源矩阵都成立的边界:ContractBlock 和 plugin skill 固定声明 frozen axes，Host 不能把 quantification 改写成当前工具能够搜索的参数。它把“有 truth 时的胜出配置只在当前蛋白矩阵切片内比较、无 truth 时不产生 winner”讲明白。

若维护者确实需要读到具体的冻结值,建议的最小形态是:上游输出 requested/effective 双层 provenance，由新的 feature-level MCP input contract 显式读取。不能从蛋白矩阵或命令默认值反推。

## 8. 结论

方向(把 quant 选择纳入 agentic)与维护者一致且有 maxlfq/directlfq 证据支撑。补测已完成,把候选范围**收敛**到这两种:piBAQ/topn/sum 在 LFQ 上从不胜出、在 TMT 上根本不适用,不进候选空间。

**run 归一化并入同一批**(第 7 节):它与 quant 共用同一个架构障碍与同一套改动,不单独立项。

剩下的仍是**真架构改动**(第 3 节的三处 + 输入契约),且成本-收益比因候选集只有两种而变窄。建议:**架构实现待维护者(Yasset)点头再开工**,本轮不实现代码。
