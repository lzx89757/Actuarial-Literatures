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

HMC 算法源于物理中研究的小球在光滑平面上的运动轨迹，通过物理模型的建模思路来进行随机模拟的模型建构

算法介绍：

* $\phi$ ：中间变量

* $\varepsilon$ ：梯度（时间）

* ：混淆矩阵

* ：参数向量（x 左边）

* ：后验密度（值）

算法步骤：

1. 首先
2. 然后
3. 其次

----

### Week 2: 随机性准备金评估模型
参考文献：[Stochastic loss reserving using bayesian MCMC models.pdf](https://github.com/lzx89757/Actuarial-Literatures/blob/master/papers/Stochastic%20loss%20reserving%20using%20bayesian%20MCMC%20models.PDF)

随机性准备金评估模型主要包含 Mack 模型、Bootstrap ODP 模型、相关链锑模型、相关增量趋势模型和结案率变化模型。

#### 数据说明

4个险种(Commercial Auto, Personal Auto, Workers Compensation, Other
Liability)中，每个险种选择 50 家保险公司的[准备金数据](http://www.casact.org/research/index.cfm?fa=loss_reserves_data)。

#### 模型设定


* **Mack model(1993, 1994)**
  Mack 模型是链锑法的推广，假设

* **Bootstrap ODP**
  假设增量赔款服从过离散的泊送分布，可以运用 GLM 进行估计，再运用 Bootstrap 抽样计算预测值的方差








* Corelated Chain-Ladder(CCL) Model - 相关链锑模型
  假设累积赔款
  服从对数正态分布

  $\mu_{w,d}=\alpha_{w} + \beta_{\alpha}+\rho(\log(C_{w,d})-\mu_{w,d})$  

* Leveled Chian Ladder(LCL) Model - 水平链

其中 $w,d$ 分别对应事故年和进展年

* CIT（相关增量趋势模型）
* CSR（结案率变化模型）


-------------------
### Week 3: 极值理论与广义帕累托分布
* [Estimating extreme tail risk measures with generalized Pareto distribution.pdf](https://github.com/lzx89757/Actuarial-Literatures/blob/master/papers/Estimating%20extreme%20tail%20risk%20measures%20with%20generalized%20Pareto%20distribution.pdf)




-------------------
### 

* [A max-stable process model for rainfall extremes at different accumulation durations 2016 Weather and Climate Extremes.pdf](https://github.com/lzx89757/Actuarial-Literatures/blob/master/papers/A%20max-stable%20process%20model%20for%20rainfall%20extremes%20at%20different%20accumulation%20durations%202016%20Weather%20and%20Climate%20Extremes.pdf)

### 
* [Tukey max-stable processes for spatial extremes.pdf](https://github.com/lzx89757/Actuarial-Literatures/blob/master/papers/Tukey%20max-stable%20processes%20for%20spatial%20extremes.pdf)

### 
* [A flexible dependence model for spatial extremes.pdf](https://github.com/lzx89757/Actuarial-Literatures/blob/master/papers/A%20flexible%20dependence%20model%20for%20spatial%20extremes.pdf)




