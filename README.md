# 风险管理与精算学术研讨会

> 时间：2017 年 2 月 - 2017 年 6 月
>
> 地点：中国人民大学明德主楼

研讨会主要内容包括：

* 贝叶斯方法

* 极值理论(GEV theory)

* 泊松过程(Poisson piont process)

* Max-stable process

* 空间模型

-------------------
### Week 1: Hamiltonian Monte Carlo 算法

参考文献：

- [Faster estimation of Bayesian models in ecology using Hamiltonian Monte Carlo.pdf](https://github.com/lzx89757/Seminar-2017/blob/master/week%201/Faster%20estimation%20of%20Bayesian%20models%20in%20ecology%20using%20Hamiltonian%20Monte%20Carlo.pdf) 及其 [R 代码](https://github.com/colemonnahan/gradmcmc)
- [贝叶斯方法在 Rstan 中的应用](https://github.com/lzx89757/Introduction-to-Rstan)
- Bayesian data analysis - Chapter 11
- [The No-U-Turn Sampler.pdf](https://github.com/lzx89757/Seminar-2017/blob/master/week%201/The%20No-U-Turn%20Sampler.pdf)
- [讲稿 PPT](https://github.com/lzx89757/Seminar-2017/blob/master/week%201/Faster%20estimation%20of%20Bayesian%20models(PPT).pptx)

内容主要涉及到了 Rejection sampling、Importance sampling、
Metropolis-Hastings 以及 Hamiltonian Monte Carlo 算法的具体步骤。Hamiltonian Monte Carlo 算法是 MCMC 算法的一种，相当于在抽样过程中增加了一个水平方向的趋势，使得 HMC 算法在抽样过程中能够遍历整个值域

----

### Week 2: 随机性准备金评估模型

#### 参考文献：

1. [Stochastic loss reserving using bayesian MCMC models.pdf](https://github.com/lzx89757/Actuarial-Literatures/blob/master/papers/Stochastic%20loss%20reserving%20using%20bayesian%20MCMC%20models.PDF)
2. A Practitioner’s Introduction to Stochastic Reserving.pdf

数据包含四个险种(Commercial Auto, Personal Auto, Workers Compensation, Other
Liability)中，每个险种选择 50 家保险公司的[流量三角形数据](http://www.casact.org/research/index.cfm?fa=loss_reserves_data)。数据的形式如下：

![editpathrtools](https://raw.githubusercontent.com/lzx89757/Seminar-2017/master/pictures/reserving%20data.png)

随机性准备金评估模型主要包含 **Mack 模型**、**Bootstrap ODP 模型**、**CCL模型**、**LCL模型**、**CIT模型**、**LIT模型**和**CSR模型**。模型之间的关系如下：

![editpathrtools](https://raw.githubusercontent.com/lzx89757/Seminar-2017/master/pictures/mcmc%20model.png)

####  一、Multiplicative Chainladder - Mack model(1993, 1994)

Mack 模型是，也称之为多元链锑法模型，是传统链锑法的推广，其 BLUE 的估计值与链锑法得到的结果等价。该方法将累积赔款 $\tilde{C}_{w,d+1}$ 作为随机变量，属于随机性准备金评估模型的一种，模型设定如下：

$$
\begin {aligned} 
  &\text{E}[\tilde{C}_{w,d+1}|C_{w,1},...,C_{w,d}]=C_{w,d}\cdot f_{d}\\
  &\text{Var}[\tilde{C}_{w,d+1}|C_{w,1},...,C_{w,d}]=C_{w,d}\cdot \alpha_{d}^{2}
\end {aligned}
$$

其中，事故年 $w=2,...,K$ 的累积赔款预测值为

$$
\begin {aligned} 
\hat{C}_{w,K}&=C_{w,k+1-w}\cdot \hat{f}_{K+1-w}\cdot \cdots \cdot \hat{f}_{K-1}\\
 \hat{f}_{d}&=\frac{\sum_{w=1}^{K-d}C_{w,d+1}}{\sum_{w=1}^{K-d}C_{w,d}}
\end {aligned}
$$

运用 R 软件的 *ChainLadder* 包可以计算得到累积赔款 $$\text{SD}[\tilde{C}_{w,K}]​$$ 和 $$\text{SD}[\sum_{w=2}^{K}\tilde{C}_{w,K}]​$$ 的标准差，也可以通过显示表达式计算模型的预测均方误差 $\text{MSEP}​$ (bootstrap 方法同样适用)。

**Mack 模型的缺陷在于：**

1. 假设同一事故年在不同进展年之间是相互独立的 - 随着时间的推移，进展年的累积赔款可能服从不同的部分

2. 参数过多，可能存在过拟合现象 - $\alpha_{d}$ 在不同的进展年下都不同
3. 模型可以处理负值或者零值，但是对于稀疏数据不适用
4. 对于长尾或者厚尾业务，方差参数估计有困难

#### 二、Over-Dispersed Poisson Model - ODP model (2002)

ODP 模型也称之为过离散泊松模型。该方法假设**增量赔款**服从过离散的泊送分布，可以运用 GLM 进行估计，再运用 Bootstrap 抽样计算预测值的标准差和均方误差。模型设定如下：
$$
\begin {aligned} 
\text{E}[\tilde{I}_{w,d}]&=\alpha_{w}\cdot \beta_{d}\\
\text{Var}[\tilde{I}_{w,d}]&=\phi\cdot \alpha_{w}\cdot\beta_{d}
\end {aligned}
$$

#### 三、The Correlated Chain-Ladder (CCL) Model - 相关链锑模型

**Bayesian Models for Incurred Loss Data:**

* The Correlated Chain-Ladder (CCL) Model
* The Leveled Chian Ladder (LCL) Model

假设累积赔款 $\tilde{C}_{w,d}$ 服从参数为 $(\mu_{w,d}, \sigma_{d})$ 的**对数正态分布**，且有 $\sigma_{1} > \sigma_{2} > \dots > \sigma_{10}$，同时假设事故年之间是累积赔款是相关的。模型设定如下：
$$
\begin {aligned} 
\mu_{1,d} &= \alpha_{1} + \beta_{d}\\
\mu_{w,d}&=\alpha_{w} + \beta_{\alpha}+\rho[\log(C_{w,d})-\mu_{w,d}]
\end {aligned}
$$
其中 $w,d$ 分别对应事故年和进展年，$\rho$ 为随机变量 $\log(\tilde{C}_{w-1,d})$ 与 $\log(\tilde{C}_{w,d})$ 的相关系数。模型的先验分布为：
$$
\begin {aligned} 
&\alpha_{w}\sim \text{normal}(\log(\text{Premium_{w}})+logelr,\sqrt{10})\\
&logelr \sim \text{uniform}(-1,0.5)\\
&\rho \sim \text{uniform}(-1,1)\\
&\beta_{d}\sim \text{uniform}(-5,5)\\
&\sigma_{d} = \sum_{i=d}^{10}\alpha_{i}\\
&\alpha_{i}\sim \text{uniform}(0,1)
\end {aligned}
$$

另外，当 $\rho = 0$ 时，**The Leveled Chian Ladder (LCL) Model - 分层链锑模型** 是 CCL 模型的特例。

#### 四、The Correlated Incremental Trend (CIT) Model - 相关增量趋势模型 

Bayesian Models for Paid Loss Data
* The Correlated Incremental Trend (CIT) Model
* The Leveled Incremental Trend (LIT) Model
* The Changing Settlement Rate (CSR) Model

假设增量赔款 $\tilde{I}_{w,d}$ 服从**偏正态分布** (偏正态分布等价与对数正态-正态分布)，且事故年之间的增量赔款具有相关性，同时伴随着日历年的趋势效应，即 **CIT 模型**设定如下：
$$
\begin {aligned} 
&\tilde{I}_{1,d}\sim \text{normal}(Z_{1,d},\delta)\\
&\tilde{I}_{w,d}\sim \text{normal}(Z_{w,d}+\rho\cdot(\tilde{I}_{w-1,d}-Z_{w-1,d})\cdot e^{\tau},\delta) \quad \text{for} \ \ \ w > 1\\
&Z_{w,d}\sim \text{lognormal}(\mu_{w,d},\sigma_{d})\quad\quad \sigma_{1}<...<\sigma_{10}\\
&\mu_{w,d}=\alpha_{w}+\beta_{d}+\tau \cdot(w+d-1)
\end {aligned}
$$
其中，$\rho$ 表示增量赔款 $\tilde{I}_{w,d}$ 和 $\tilde{I}_{w-1,d}$ 的相关系数。模型的先验分布为：
$$
\begin {aligned} 
&\sigma_{1}^{2}\sim \text{uniform}(0,0.5)\\
&\sigma_{d}^{2}\sim \text{uniform}(\sigma_{d-1}^{2},\sigma_{d-1}^{2} + 0.1)\\
&\alpha_{w} \sim \text{normal}(\log(\text{Premium_{w}})+logelr,\sqrt{10})\\
&logelr \sim \text{uniform}(-1,0.5)\\
&\rho \sim \text{uniform}(-1,1)\\
&\beta_{d} \sim \text{uniform}(0,10)\quad \text{for}\quad d=1,...,4\\
&\beta_{d} \sim \text{uniform}(0,\beta_{d-1})\quad \text{for}\quad d>4\\
&\tau \sim \text{normal}(0,0.0316)\\
&\delta \sim \text{uniform}(0,\text{Average Premium})
\end {aligned}
$$
其中，当 $\rho=0$ 时，**LIT 模型**是 **CIT 模型**的特例。

#### 五、 The Changing Settlement Rate (CSR) Model - 结案率变化模型

假设累积赔款 $\tilde{C}_{w,d}$ 服从参数为 $(\mu_{w,d}, \sigma_{d})$ 的**对数正态分布**，且有 $\sigma_{1} > \sigma_{2} > \dots > \sigma_{10}$，同时考虑案件处理速度的逐年变化。模型设定如下：
$$
\mu_{w,d}=\alpha_{w} + \beta_{\alpha}\cdot(1-\gamma)^{w-1}
$$
其中 $w,d$ 分别对应事故年和进展年。模型的先验分布为：
$$
\begin {aligned} 
&\alpha_{w} \sim \text{normal}(\log(\text{Premium}_{w})+logelr,\sqrt{10})\\
&logelr \sim \text{uniform}(-1,0.5)\\
&\gamma \sim \text{normal}(0,0.025)\\
&\beta_{d} \sim \text{uniform}(-5,5) \quad \beta_{10}=0 \\
&\sigma_{d} = \sum_{i=d}^{10}\alpha_{i}\\
&\alpha_{i}\sim \text{uniform}(0,1)
\end {aligned}
$$

#### 六、结论

* For Incurred Loss Data：CCL＞LCL＞Mack
* For Paid Loss Data：CSR＞CIT≈LIT＞Bootstrap ODP＞Mack

-------------------
### Week 3: 极值理论与广义帕累托分布
参考文献：
* [Estimating extreme tail risk measures with generalized Pareto distribution.pdf](https://github.com/lzx89757/Actuarial-Literatures/blob/master/papers/Estimating%20extreme%20tail%20risk%20measures%20with%20generalized%20Pareto%20distribution.pdf)

文章介绍了极值尾部风险的度量，主要分为两部分：
* 提出一种新的广义帕累托分布 (GPD) 的估计方法
* 用于得到风险度量值
* 可以运用 POT package 或者 optim 等非线性规划的函数


-------------------
### 

* [A max-stable process model for rainfall extremes at different accumulation durations 2016 Weather and Climate Extremes.pdf](https://github.com/lzx89757/Actuarial-Literatures/blob/master/papers/A%20max-stable%20process%20model%20for%20rainfall%20extremes%20at%20different%20accumulation%20durations%202016%20Weather%20and%20Climate%20Extremes.pdf)

### 
* [Tukey max-stable processes for spatial extremes.pdf](https://github.com/lzx89757/Actuarial-Literatures/blob/master/papers/Tukey%20max-stable%20processes%20for%20spatial%20extremes.pdf)

### 
* [A flexible dependence model for spatial extremes.pdf](https://github.com/lzx89757/Actuarial-Literatures/blob/master/papers/A%20flexible%20dependence%20model%20for%20spatial%20extremes.pdf)




