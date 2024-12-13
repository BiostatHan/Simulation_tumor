# 肿瘤临床试验模拟与分析工具包

本项目提供了两个主要函数 `simulation` 和 `analyze`，用于模拟肿瘤临床试验数据并进行相关分析。以下是详细的函数介绍、参数说明及内部实现逻辑。

---

## Simulation 函数

### 1.1 函数介绍

`simulation` 函数用于模拟肿瘤临床试验的数据。通过设定试验和对照组的客观响应率（ORR）、显著性水平、检验效能等参数，生成一个包含受试者基本信息、响应状态及时间指标的数据框。生成的数据框包含受试者的编号、分组、入组时间、是否达到响应、初始响应状态、达到响应的时间、响应持续时间、无进展生存时间以及不同访视时间点的疾病状态。

#### 1.1.1 传入参数

- `p1`：试验组ORR，默认值为0.5。
- `p2`：对照组ORR，默认值为0.4。
- `alpha`：显著性水平，默认值为0.05。
- `power`：检验效能，默认值为0.8。
- `dropout`：预期脱落率，默认值为0.2。
- `enrollment_period`：入组期（月），默认值为12。
- `followup_period`：随访期（月），默认值为24。
- `assessment_interval`：评估间隔（月），默认值为6。
- `cr_pr_ratio`：CR/PR比例，默认值为 `c(0.3, 0.7)`。
- `cr_pfs`：CR的无进展生存时间均值，默认值为15。
- `pr_pfs`：PR的无进展生存时间均值，默认值为10。
- `non_responder_pfs`：非响应者的无进展生存时间均值，默认值为5。
- `ttr`：达到响应的时间均值，默认值为2。
- `seed`：随机种子，默认值为12345。

#### 1.1.2 返回数据

函数返回一个数据框，包含以下主要字段：

- `SubjID`：受试者编号。
- `Group`：分组（Treatment/Control）。
- `InitialResponse`：初始响应状态（CR/PR/PD）。
- `EnrollmentTime`：入组时间。
- `TTR`：达到响应的时间（Time to Response）。
- `DoR`：响应持续时间（Duration of Response）。
- `PFS`：无进展生存时间（Progression-Free Survival）。
- `Month_X`：不同随访时间点的疾病状态（例如 `Month_0`, `Month_6`, `Month_12`, 等）。

### 1.2 数据生成逻辑

#### 1.2.1 样本量计算

基于以下参数计算每组所需的样本量：

- `p1` 和 `p2`：试验组和对照组的预期ORR。
- `alpha`：显著性水平。
- `power`：检验效能。
- `dropout`：预期脱落率。

使用标准统计公式计算每组的样本量，并考虑脱落率调整总体样本量。

#### 1.2.2 入组时间生成

在设定的入组期内，均匀分布随机生成每位受试者的入组时间：

\[
\text{EnrollmentTime} \sim \text{Uniform}(0, \text{enrollment_period})
\]

#### 1.2.3 响应状态生成

根据预设的响应率 `p1` 和 `p2`，分别为试验组和对照组生成二项分布的响应状态：

- **试验组**：`InitialResponse` ~ Binomial(1, `p1`)。
- **对照组**：`InitialResponse` ~ Binomial(1, `p2`)。

对于响应的受试者，进一步根据 `cr_pr_ratio` 随机分配为完全缓解（CR）或部分缓解（PR）状态；未响应的标记为进展疾病（PD）。

### 1.3 时间指标生成

#### 1.3.1 响应者

对于达到响应的受试者（CR或PR），生成以下时间指标：

- **TTR（Time to Response）**：服从指数分布，中位数约为2个月。
  
  \[
  \text{TTR} \sim \text{Exponential}(\text{rate}=1/2)
  \]

- **PFS（Progression-Free Survival）**：
  - **CR患者**：

    \[
    \text{PFS} \sim \text{Exponential}(\text{rate}=1/\text{cr_pfs})
    \]

  - **PR患者**：

    \[
    \text{PFS} \sim \text{Exponential}(\text{rate}=1/\text{pr_pfs})
    \]

- **DoR（Duration of Response）**：计算为

  \[
  \text{DoR} = \text{PFS} - \text{TTR}
  \]

#### 1.3.2 非响应者

对于未达到响应的受试者（PD），生成以下时间指标：

- **PFS（Progression-Free Survival）**：服从指数分布，中位数约为5个月。

  \[
  \text{PFS} \sim \text{Exponential}(\text{rate}=1/5)
  \]

### 1.4 疾病状态生成

在每个随访时间点，根据以下规则确定受试者的疾病状态：

#### 1.4.1 对于响应者

- 如果当前时间 < `EnrollmentTime`：标记为 `"Not Enrolled"`。
- 如果 `EnrollmentTime` <= 当前时间 < (`EnrollmentTime` + `TTR`)：标记为 `"SD"`（稳定疾病）。
- 如果 (`EnrollmentTime` + `TTR`) <= 当前时间 <= (`EnrollmentTime` + `PFS`)：标记为初始响应状态（`"CR"` 或 `"PR"`）。
- 如果当前时间 > (`EnrollmentTime` + `PFS`)：标记为 `"Exited/Died"`。

#### 1.4.2 对于非响应者

- 如果当前时间 < `EnrollmentTime`：标记为 `"Not Enrolled"`。
- 如果 `EnrollmentTime` <= 当前时间 <= (`EnrollmentTime` + `PFS`)：标记为 `"PD"`（进展疾病）。
- 如果当前时间 > (`EnrollmentTime` + `PFS`)：标记为 `"Exited/Died"`。

## Analyze 函数

### 2.1 函数介绍

`analyze` 函数用于分析 `simulation` 函数生成的临床试验数据。在指定的随访时间点，分析各组的疗效指标，包括客观响应率（ORR）、达到响应时间（TTR）、响应持续时间（DoR）以及无进展生存时间（PFS）。函数还会进行统计检验，如Fisher精确检验和生存分析。

#### 2.1.1 传入参数

- `data`：包含临床试验数据的数据框，由 `simulation` 函数生成。
- `timepoint`：指定的随访时间点（以月为单位），用于分析该时间点的疾病状态。

#### 2.1.2 返回数据

函数返回一个包含以下内容的列表：

- `status_counts`：各组在指定时间点的疾病状态人数统计。
- `orr`：试验组和对照组的客观响应率（ORR）。
- `fisher_p`：ORR差异的Fisher精确检验p值。
- `dor`：DoR分析结果，包括中位数、生存时间置信区间、风险比（HR）及其置信区间、p值。
- `ttr`：TTR分析结果，包括中位数、生存时间置信区间、风险比（HR）及其置信区间、p值。
- `pfs`：PFS分析结果，包括中位数、生存时间置信区间、风险比（HR）及其置信区间、p值。

### 2.2 计算细节

#### 2.2.1 获取时间点

- 使用 `paste0` 函数构建包含时间点的列名（例如，`Month_12`）。
- 检查该列名是否存在于数据框中，如果不存在，则停止执行并返回错误信息。

#### 2.2.2 ORR 计算

- 从数据框中提取指定时间点的疾病状态。
- 使用 `group_by` 和 `summarise` 计算试验组和对照组达到完全缓解（CR）或部分缓解（PR）的患者比例。
- 使用 `fisher.test` 进行Fisher精确检验，以评估两组响应率的差异是否具有统计学意义。

#### 2.2.3 生存分析

进行三项生存分析：

1. **DoR（响应持续时间）分析**：
   - 筛选有 `DoR` 值的患者。
   - 使用K-M曲线（`survfit`）和Cox回归（`coxph`）计算中位DoR、生存时间置信区间、风险比（HR）及其置信区间。
   - 使用Log-Rank检验（`survdiff`）计算p值。

2. **TTR（达到响应时间）分析**：
   - 筛选有 `TTR` 值的患者。
   - 同样使用K-M曲线、Cox回归和Log-Rank检验计算相关指标。

3. **PFS（无进展生存时间）分析**：
   - 使用所有已入组患者的 `PFS` 数据进行分析。
   - 计算中位PFS、生存时间置信区间、风险比（HR）及其置信区间。
   - 使用Log-Rank检验计算p值。

所有统计结果将通过 `pander` 包进行格式化输出。

## 示例

以下是如何使用 `simulation` 和 `analyze` 函数的示例：

```R
# 加载必要的包
library(tidyverse)
library(survival)
library(pander)

# 模拟数据
set.seed(12345)
sim_data <- simulation(
  p1 = 0.5,
  p2 = 0.4,
  alpha = 0.05,
  power = 0.8,
  dropout = 0.2,
  enrollment_period = 12,
  followup_period = 24,
  assessment_interval = 6,
  cr_pr_ratio = c(0.3, 0.7),
  cr_pfs = 15,
  pr_pfs = 10,
  non_responder_pfs = 5,
  ttr = 2,
  seed = 12345
)

# 分析指定时间点（例如12个月）的数据
analysis_results <- analyze(data = sim_data, timepoint = 12)
