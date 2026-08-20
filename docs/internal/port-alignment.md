# R → Python 移植对齐基准

> 状态:**已执行并归档**(2026-07-16)。本文档记录 mokume 中四个「参考实现在 R、内核重写在 Python」的组件与原 R 包的数值核对结果,数据、脚本、命令与结论均可复现。

移植件采用与否的判据不是「像」,而是:**逐行移植的必须数值一致,非逐行移植的必须如实标注、只声称同族/同量级/同方向**。本文档对四项分别给出裁决。

## 0. 复现环境与数据

### 0.1 环境

临时 conda env(带 R + rpy2),本仓库运行时**不依赖**它,它只为本对齐存在:

```bash
source /home/shenyufei/miniforge3/etc/profile.d/conda.sh
conda activate mokume_portalign2
export PYTHONPATH=/home/shenyufei/Git-repository/Bigbio/mokume/python
```

| 组件 | 版本 |
| --- | --- |
| R | 4.5.3 (2026-03-11) |
| `qvalue` | 2.42.0 |
| `fdrtool` | 1.2.18 |
| `locfdr` | 1.1.8 |
| `samr` | 3.0.1 |
| Python | 3.12 |
| numpy / scikit-learn | 2.4.6 / 1.9.0 |

`permFDP` 未在本环境安装,也未在 CRAN / Bioconductor 找到该名称的包 —— **不可得**,第 4 项因此只与 `samr` 对齐。

### 0.2 数据

全部为**真实 spike-in 数据**,不是模拟:mokume 网格已产出的 DE 表
`/data/shenyufei/Bigbio_data/PXD_spike_in/PXD*/mokume/grid_results/de/*.tsv`
(列 `ProteinName, log2FC, pvalue, adj_pvalue, ...`)。

抽样规则:13 个 PXD,每个按文件名排序后**等间距取 6 张**(确定性,无随机),第 1/2 项再要求 ≥100 个有限 p 值,第 3 项要求 ≥200 个有限 log2FC。

| 覆盖 | 值 |
| --- | --- |
| 数据集 | 13 个(PXD000279 / 003881 / 007683_LFQ / 015261 / 020815 / 026600 / 028735 / 040449 / 046444 / 062685 / 066672 / 070049 / 070151) |
| DE 表 | 77 张(第 1/2 项),71 张(第 3 项) |
| 每张表假设数 | 165 ~ 14731 |
| R 侧 π0 跨度 | 0.0056 ~ 1.0000 |

### 0.3 脚本与原始输出

脚本与逐表原始数值落盘在 scratchpad(**不入库**,因其依赖上面的临时 env 与本机数据路径):

```text
/tmp/claude-1001/-home-shenyufei-Git-repository-Bigbio/a3fa24d6-6dc3-419e-bd05-82ff80cff97c/scratchpad/
  align_common.py                # 取数helper
  align_01_qvalue_pi0.py   -> align_01_qvalue_pi0.tsv   / .log   # 第 1、2 项
  align_02_effect_gate.py  -> align_02_effect_gate.tsv  / .log   # 第 3 项
  align_03_samr.py         -> align_03_samr.tsv         / .log   # 第 4 项
  pxd028735_setup.py             # PXD028735 夹具(第 4 项与置换验证共用)
```

---

## 1. Storey q-value vs R `qvalue::qvalue` —— **对齐,采用**

**对象**:`mokume.analysis.adaptive_fdr.qvalues`。

**协议**:分两层比,才能把 π0 估计的误差和 q 值计算本身的误差分开。

1. **把 π0 固定成 R 的估计值**,两边只跑 q 值那一步 —— 这是对 `qvalues()` 本身的干净检验。
2. **端到端**(各自估各自的 π0),这是使用者实际拿到的东西。

**结果**(72 张可比表):

| 比较 | max \|diff\| | median \|diff\| |
| --- | --- | --- |
| q 值(π0 固定为 R 的) | 2.220e-16 | 5.551e-17 |
| q 值(端到端) | 7.155e-05 | 4.748e-06 |
| q ≤ 0.05 的发现数差异 | **0**(72/72 张完全相同) | 0 |

**结论**:π0 固定后差异在 **机器精度**(1e-16 量级 = double 舍入)—— `qvalues()` 是 `qvalue::qvalue` 的**逐行等价移植**。端到端 1e-5 量级的残差**全部来自 π0**(见第 2 项),且不改变任何一张表在 q ≤ 0.05 下的发现集合。**采用**。

## 2. π0 smoother vs R `qvalue::pi0est` —— **对齐,采用**

**对象**:`mokume.analysis.adaptive_fdr.estimate_pi0(method="smoother")`,对齐时 **`conservative_bound=False`**(下界是 mokume 自己加的护栏,R 没有对应物,开着比就不是同一个估计量了)。

R 侧:`pi0est(p, pi0.method="smoother")$pi0`。

**结果**(72 张可比表):

| 指标 | 值 |
| --- | --- |
| π0 \|diff\| max | 7.247e-05 |
| π0 \|diff\| median | 4.772e-06 |

**结论**:**采用**。残差 1e-5 量级,来源是两边把 `smooth.spline(..., df=3)` 的有效自由度钉到 3 的方式不同 —— R 用自己的 Reinsch 实现,mokume 用对 log(λ) 的 60 次二分。这是**收敛容差**,不是模型差异;它比 π0 本身小 4~5 个数量级,且如第 1 项所示对发现集合零影响。

### 2.1 R 拒绝、移植件不拒绝的 5 张表(必须记录)

77 张里有 5 张 R **直接报错**,无法比较:

| 数据集 | R 报错 | 移植件返回的 π0(raw) | 下界护栏 |
| --- | --- | --- | --- |
| PXD020815 | `pi0est`: estimated pi0 <= 0 | 2.2e-308(≈0,被 clamp 到 float 最小正数) | **救回 0.0061** |
| PXD046444 | `pi0est`: estimated pi0 <= 0 | 2.2e-308 | **救回 0.0303** |
| PXD062685 | `pi0est`: estimated pi0 <= 0 | 2.2e-308 | **救回 0.0696** |
| PXD000279 | `smooth.spline`: 不允许缺失/无限值 | 2.2e-308 | 未生效(曲线本身触 0,退化表) |
| PXD015261 | `smooth.spline`: 不允许缺失/无限值 | 0.1372 | 未生效 |

**这三行是 `conservative_bound` 迄今最强的真实数据依据**:smoother 把样条外推到 λ→1 时确实会落到 0 以下,R 的做法是**报错罢工**,移植件的做法是 clamp 到 float 最小正数(等于给自适应 FDR 一个无限大的拒绝预算,比 R 更危险),而下界护栏把它抬回 π0(λ) 曲线支持得住的值。

### 2.2 关于 `conservative_bound` 的诚实修正

`adaptive_fdr.pi0_lower_bound` 的 docstring 写:「on simulated two-group p-values it binds on roughly 1% of datasets and then raises pi0 by well under 1%」。**这句话在真实 spike-in 数据上不成立**:

| 指标 | docstring 的模拟数据说法 | 本次真实数据实测(77 张) |
| --- | --- | --- |
| 生效频率 | ~1% | **31%**(24/77) |
| 生效时抬升 π0 | 「well under 1%」 | 绝对值 median +0.0098 / **max +0.0902**;相对 median ×1.054 / **max ×3.258** |

方向没错(护栏只会让 FDR 更保守,不会更激进),量级差了一个数量级以上:真实蛋白组 DE 的 p 值分布比模拟的两组正态更不规整,smoother 外推离开曲线的情况远比模拟常见。源码 docstring 已更新为本节的真实数据结果。

---

## 3. 效应量门 vs R `fdrtool` / `locfdr` —— **同族/同方向,采用;但不是逐行移植**

### 3.1 先说清楚它不是什么

`mokume.analysis.effect_size_gate.estimate_effect_size_gate(method="mixture")` **不是 `fdrtool` 的移植**,一行都不是。

| | mokume | fdrtool | locfdr |
| --- | --- | --- | --- |
| 零分布 | 对 \|log2FC − median\| 拟合**两成分高斯混合**,零成分即其中较窄的那个 | 自己的半正态/经验零分布拟合 | 自己的经验零(`nulltype`,默认 ML) |
| 门的定义 | 两成分后验概率交点(该混合模型下的 local-fdr = 0.5 边界) | 其 lfdr 曲线 | 其 fdr 曲线 |
| 实现来源 | sklearn `GaussianMixture` | Strimmer 的 grenander/半正态内核 | Efron 的多项式 logit 拟合 |

三者**同属 local-fdr 家族**,但零分布模型不同,**所以不应期待数值相等**。本节能问的只有三件事:同族?同量级?同方向?

### 3.2 对齐协议

三者统一放到 **log2FC 尺度**、统一取 **local-fdr = 0.5** 边界:

- mokume:`estimate_effect_size_gate(x, method="mixture")`
- fdrtool:`fdrtool(x - median(x), statistic="normal")$lfdr` → 门 = lfdr ≤ 0.5 的最小 \|x − median(x)\|
- locfdr:`locfdr(z)$fdr`,z = (x − median)/(1.4826·MAD) → 门 = fdr ≤ 0.5 的最小 \|z\|,**再乘回 1.4826·MAD** 回到 log2FC 尺度(线性映射,不改变排序)

### 3.3 结果(71 张表)

门的分布:

| 估计量 | 可用表数 | median | IQR | range |
| --- | --- | --- | --- | --- |
| mokume mixture | **71/71** | 0.443 | [0.221, 0.873] | [0.009, 5.299] |
| fdrtool lfdr50 | 66/71 | 0.996 | [0.562, 2.282] | [0.005, 31.367] |
| locfdr lfdr50 | **50/71** | 0.510 | [0.307, 1.175] | [0.061, 4.574] |

两两一致性:

| 配对 | median 比值 | p10–p90 | 落在 2 倍以内 | Spearman ρ |
| --- | --- | --- | --- | --- |
| mokume / **locfdr** | **1.019** | [0.292, 1.815] | **82%** | 0.765 (p=9.5e-11, n=50) |
| mokume / fdrtool | 0.362 | [0.180, 1.362] | 29% | 0.637 (p=9e-09, n=66) |
| **fdrtool / locfdr**(两个 R 参考互比) | **2.631** | [1.206, 6.079] | **31%** | 0.856 (p=9e-15, n=48) |

### 3.4 裁决

- **同族**:是。三者都能在同一批 log2FC 上定出 local-fdr = 0.5 边界,且 mokume 从未退化到 `fallback=0.5`(0/71 次触发回退,即 71 张门全部是数据推出来的)。
- **同方向**:是。ρ = 0.64(vs fdrtool)/ 0.77(vs locfdr),均高度显著。
- **同量级**:对 locfdr 是(中位比值 1.02,82% 落在 2 倍以内);对 fdrtool 系统性偏松 ~2.8 倍(mokume 门更小,53/66 张更低)。
- **关键背景**:**两个 R 参考自己就差 2.6 倍**(fdrtool / locfdr 中位比值 2.631,只有 31% 落在 2 倍以内)。也就是说,mokume 与 locfdr 的分歧(1.02)**小于两个 R 参考之间的分歧**(2.63)。「门定在哪」这件事本身在 R 生态里就没有唯一答案,mokume 落在 R-vs-R 的散布包络内。
- **鲁棒性**:locfdr 在 21/71 张(30%)上直接失败(`CM estimation failed, middle of histogram non-normal`、`solve(G0)` 奇异、`lm.fit` 报错),fdrtool 在 5/71 张(7%)上失败;mokume 0 张失败。

**结论:采用**,并如实标注为「同族的独立实现,不是 fdrtool/locfdr 的逐行移植,数值不与任一 R 参考精确相等,亦不应精确相等」。`effect_size_gate.py` 模块 docstring 已经把这点写对了(「in the spirit of the R `fdrtool` / `locfdr` local-fdr thresholding, but native and dependency-free」「methodological lineage, not a direct citation」),本节为其补上真实数据证据。

---

## 4. 汇总

| # | 移植件 | R 参考 | 是否逐行移植 | 数值差异 | 裁决 |
| --- | --- | --- | --- | --- | --- |
| 1 | `adaptive_fdr.qvalues` | `qvalue::qvalue` | **是** | π0 固定时 max 2.2e-16(机器精度);q≤0.05 发现集 72/72 完全相同 | **采用** |
| 2 | `adaptive_fdr.estimate_pi0(smoother)` | `qvalue::pi0est` | **是** | max 7.2e-05 / median 4.8e-06(样条 df 二分容差) | **采用** |
| 3 | `effect_size_gate.estimate_effect_size_gate` | `fdrtool` / `locfdr` | **否**(同族独立实现) | vs locfdr 中位比 1.02、ρ=0.77;vs fdrtool 中位比 0.36、ρ=0.64;**两个 R 参考互比 2.63** | **采用**(同族/同量级/同方向;不声称数值相等) |
| — | `permFDP` 对齐 | — | — | — | **不可得**(该包在本环境与 CRAN/Bioconductor 均未找到) |
