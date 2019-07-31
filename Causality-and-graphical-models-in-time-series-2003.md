## Causality and graphical models in time series （2003）

### Introduction 

​		过去几年里，人们对于图模型，尤其是基于有向无环图（directed acyclic graphs, DAG）作为描述和推断因果关系的一个基础框架越来越感兴趣。这种新的图形化方法与形式化因果关系的概念有关，如Neyman和Rubin的潜在反应模型以及路径分析或结构方程模型。后一概念（路径分析、结构方程模型）被经济学家用来描述随时间变化的系统均衡分布。

​		本文将讨论结构方程模型的动态解释背后的思想，并讨论**利用时间来进行因果推断的图模型**。注意本文使用的因果关系不是基于干预的，而是基于一个明显的事实：**在时间上，结果无法超越其原因**。若变量在不同时间点被观测到，则可以检验各变量在不同时滞下的（偏）相关性来推断此系统中变量间是否相关。用同样的方法，可以在图模型中确定有向边的方向，从而假设先验的变量顺序，如多元数据的DAG。由于我们对变量在时间上多次测量，故这是一个时间序列。注意本文仅关注一个长期平稳的多元时间序列，而不是面板数据）。

​		因果关系的概念本文使用Granger因果关系的概念，用自然时间顺序表示因果顺序，即**若A为B的Granger原因，则用含A的信息去预测B比不用更好**。具体来说，本文将按如下顺序展开：

- **Section 2** : Eichler 曾用Granger因果关系的定义去定义时间序列的因果关系图，这一部分将这些图和时间序列链图（time series chain graphs）、时间序列偏相关图（partial correlation graphs）放在一起讨论

- **Section 3** : 讨论马尔可夫性（Markov properties）

- **Section 4** : 讨论这些图的统计推断

- **Section 5** : 两个案例分析

- **Section 6** : 结论。讨论用现有的图进行因果推断以及因果效应错误识别的可能来源

  

### Section 2 : Graphical models for multivariate time series

​		令$X=\left\{X_{a}(t), t \in \mathbb{Z}, a=1, \ldots, d\right\}$为 $d$ 元平稳序列。本文假设$X$有正定谱矩阵（positive definite spectral matrix）$f(\lambda)$，特征值有界且对于$\lambda \in [-\pi,\pi]$其一致地不等于0。又令$V=\{1,...,d\}$为指标集（set of indices）。对于$\forall A\subset V$定义$X_A=\{X_A(t)\}$为$A$中指标给出的多元子过程（subprocess）。进一步令$\overline{X}_{A}(t)=\left\{X_{A}(s), s<t\right\}$为$X_A$在 $t$ 时的历史。

​		定义 $X$ 的图模型有很多方式。这里区分两种典型的图模型。

- **第一种**：变量$X_a(t)$在 $t$ 时刻被图中单独的顶点（vertex）表示，得到经典图模型（如[定义2.1](#df2.1)引入的时间序列链图）的一般化。
- **第二种**：顶点集仅由序列的 $X_a$ 部分组成，得到了一个序列依赖结构的粗模型（coarser modelling）。这一方式会得到混合图（mixed graph），其中有向边反应Granger因果关系，而无向边代表同时期的依赖结构（dependence structure）。

这些图统称为**Granger因果关系图**（[定义2.4](#df2.4)）。此外本文还有时间序列的偏相关图，即把经典的浓度图推广到时间序列上。

​		本文只讨论线性关联和线性Granger非因果关系，因为这些线性关系可以在Hilbert随机变量空间中闭线性子空间的条件正交性来表示（看不懂没关系，看懂后面的就行）。对于随机向量 $X ,Y,Z $ ，若 $X$ 和 $Y$ 在 $ Z$ 的线性影响被移除的情况下是不相关的，则称给定 $Z$ 时 $X$ 和 $Y$ 是条件正交的，记作$X \perp Y | Z$。

> 显然这种线性方法无法得到具有部分非线性特征过程的全部因果关系，注意是统计推断的问题导致了线性限制的存在。例如，存在估计线性依赖的非参数方法，而在一般情况下，考虑到维度灾难（curse of dimensionality），非参数推断看起来难以实现。即使使用参数方法，由于大量的参数导致非线性模型往往不切实际，所以很多应用被限制在线性模型下。

作者注意到图模型方法可以用条件独立性代替条件正交性，产生强Granger因果关系，进而被扩展到非线性情况。对于高斯过程来说，图的两种意义是相同的。本章大部分结果也适用于这些一般的图，具体细节将在Section6中讨论。

#### 2.1 Time series chain graphs

​		定义图时间序列模型的第一个方法很自然地得到链图。前人介绍了基于链图经典的LWF马尔可夫性质的动态交互模型。这里我们讨论另一种方法，其根据AMP马尔可夫性质去定义图。（有关这俩性质可见另一篇文章）

> 这样定义的其中一个原因是，在时间序列链图中使用AMP马尔可夫性质和大量时间序列模型的递推结构之间存在着密切关系，这些模型可以立刻绘制出图的特征。另一个原因是AMP马尔科尔夫性质和Granger因果关系有关，而后者被用来定义Granger因果关系图。尤其是AMP马尔可夫性质允许通过对时间序列链图进行简单的聚合而得到Granger因果关系图。

<span id="df2.1">**定义2.1**</span>：【时间序列链图】平稳序列 $X$ 的时间序列链图（TSC-graph）为链图$G_{TS}=(V_{TS},E_{TS})$，其中$V_{\mathrm{TS}}=V \times \mathbb{Z}$ 且边集合（edge set）$E_{TS}$ 满足
$$
(a, t-u) \longrightarrow(b, t) \notin E_{\mathrm{TS}} \Leftrightarrow u \leq 0 \text { or } X_{a}(t-u) \perp X_{b}(t) | \overline{X}_{V}(t) \backslash\left\{X_{a}(t-u)\right\} 
$$

$$
(a, t-u)-(b, t) \notin E_{\mathrm{TS}} \Leftrightarrow u \neq 0 \text { or } X_{a}(t) \perp X_{b}(t) | \overline{X}_{V}(t) \cup\left\{X_{V \backslash\{a, b\}}(t)\right\}
$$

> 通俗解释版本：某时刻的 $a$ 不是 $t$ 时刻 $b$ 的原因，等价于该时刻大于等于 $t$ ，或给定除了该时刻的 $a$ 之外所有点的过去的情况下这俩正交；某时刻的 $a$ 和是 $t$ 时刻 $b$ 不是同时期依赖，等价于该时刻不为t，或给定所有点的过去以及 $t$ 时除这俩点以外的其他点的情况下，这俩正交

由于序列是平稳的，我们有$(a, t)-(b, t) \notin E_{\mathrm{TS}}$当且仅当$\forall s \in \mathbb{Z},(a, s)-(b,s) \notin E_{\mathrm{TS}}$。同样的时移不变性（shift invariance）同样适用于有向边。进一步注意到，上述条件保证了该过程满足 $G_{TS}$ 的成对AMP马尔可夫性质。

**例子2.2**：【向量自回归过程】假设 $X$ 是一个线性向量自回归（VAR）过程
$$
X(t)=A(1) X(t-1)+\ldots+A(p) X(t-p)+\varepsilon(t)
$$
这里误差 $\varepsilon(t)$ 独立同分布，均值为0，协方差阵为 $\Sigma$ ，若 $G_{TS}$ 表示TSC图，则其可被表示为
$$
(a, t-u) \longrightarrow(b, t) \in E_{\mathrm{TS}} \Leftrightarrow u \in\{1, \ldots, p\} \text { and } A_{b a}(u) \neq 0
$$

> $t-u$ 时刻的 $a$ 是 $t$ 时刻的 $b$ 的原因，等价于 $A_{ba}(u)$ 不为0，注意角标是 $ba$ ，说明是列对行的影响

即图中的有向边反映了时间序列的递归结构。进一步，无向边需要误差的协方差选择模型，即令 $K=\Sigma^{-1}$， 有
$$
(a, t)-(b, t) \in E_{\mathrm{TS}} \Leftrightarrow \varepsilon_{a}(t) \perp \varepsilon_{b}(t) | \varepsilon_{V \backslash\{a, b\}}(t) \Leftrightarrow K_{a b} \neq 0
$$

> $t$ 时刻 $a$ 和 $b$ 同期依赖，等价于给定 $t$ 时刻除这俩点以外所有点的误差时，这俩点的误差正交，也等价于 $K_{ab}$ 不为0

作为一个例子，可以考虑对于5元VAR(2)的过程，参数如下：
$$
A(1)=\left(\begin{array}{ccccc}{\frac{3}{5}} & {0} & {\frac{1}{5}} & {0} & {0} \\ {0} & {\frac{3}{5}} & {0} & {-\frac{1}{5}} & {0} \\ {\frac{2}{5}} & {\frac{3}{3}} & {\frac{3}{5}} & {0} & {0} \\ {0} & {0} & {0} & {-\frac{1}{2}} & {\frac{1}{2}} \\ {0} & {0} & {\frac{1}{5}} & {0} & {\frac{2}{5}}\end{array}\right), A(2)=\left(\begin{array}{ccccc}{0} & {0} & {-\frac{1}{5}} & {0} & {0} \\ {0} & {0} & {0} & {0} & {0} \\ {0} & {0} & {0} & {0} & {0} \\ {0} & {0} & {0} & {0} & {\frac{1}{3}} \\ {0} & {0} & {\frac{1}{5}} & {0} & {\frac{1}{3}} \\ {0} & {0} & {0} & {0} & {-\frac{1}{5}}\end{array}\right), \Sigma=\left(\begin{array}{ccccc}{1} & {\frac{1}{2}} & {\frac{1}{3}} & {0} & {0} \\ {\frac{1}{2}} & {1} & {-\frac{1}{3}} & {0} & {0} \\ {\frac{1}{3}} & {-\frac{1}{3}} & {1} & {0} & {0} \\ {0} & {0} & {0} & {1} & {0} \\ {0} & {0} & {0} & {0} & {1}\end{array}\right)
$$
$E_{TS}$ 中边的上述条件可以得到下图：

![fig1](https://raw.githubusercontent.com/rucwyf/Graphical-Model/master/pictures/Causality%20and%20graphical%20models%20in%20time%20series%20analysis/CAGMITSA1.png)

#### 2.2 Granger causality graphs

​		现定义贯穿全文所使用的非因果关系的定义，以及Granger因果关系图。这些混合图的顶点集仅由序列的分量组成。对于有向边，本文使用Granger因果关系的概念，对于无向边，本文使用和TSC图一样的定义。【还有很多关于因果关系的定义参见Geweke（1984），Aigner & Zellner（1988）】

<span id="df2.3">**定义2.3**</span>：【非因果关系】若满足
$$
X_{b}(t) \perp \overline{X}_{a}(t) | \overline{X}_{V \backslash\{a\}}(t)
$$
则相对于过程 $X_V$ ，$X_a$ 与 $X_b$ 无关，记作$X_{a} \nrightarrow X_{b}\left[X_{V}\right]$。进一步，若
$$
X_{a}(t) \perp X_{b}(t) | \overline{X}(t), X_{V \backslash\{a, b\}}(t)
$$
则相对于过程 $X_V$ ，$X_a$ 与 $X_b$ 偏同期不相关（partially contemporaneously uncorrelated），记作
$$
X_{a} \nsim X_{b}\left[X_{V}\right]
$$

> 通俗解释版本：给定 $X_V$ ，$X_a$ 与 $X_b$ 无关就是说给定除了 $a$ 点以外所有点 $t$ 时刻之前的信息时，$t$ 时的 $b$ 和 $t$ 时之前的 $a$ 正交；给定 $X_V$ ，$X_a$ 与 $X_b$ 偏同期不相关就是说给定 $t$ 时刻之前所有信息，以及 $t$ 时除了这俩点别的信息时，这俩点正交

注意**定义2.3**对于多元的子过程 $X_A$ 和 $X_B$ 来说也成立。在多元环境下，存在多种因果关系，如直接、间接、反馈或虚假因果关系。这些关系都可以用图来描述。

<span id="df2.4">**定义2.4**</span>：【因果关系图】平稳序列 $X$ 的Granger因果关系图是混合图 $G_C=(V,E_C)$，满足对所有 $a,b\in V$ 且 $a\ne b$ 有
$$
\begin{array}{l}{(i)\qquad a \rightarrow b \notin E_{\mathrm{C}} \Leftrightarrow X_{a} \nrightarrow X_{b}\left[X_{V}\right]} \\(ii)\qquad {a-b \notin E_{\mathrm{C}} \Leftrightarrow X_{a} \nsim X_{b}\left[X_{V}\right]}\end{array}
$$

> 通俗解释版本：在因果关系图中，$a$ 不是 $b$ 的原因等价于给定 $X_V$ ，$X_a$ 与 $X_b$ 无关；$a$ 和 $b$ 不同期相关等价于给定 $X_V$ ，$X_a$ 与 $X_b$ 偏同期不相关

简单起见，将Granger因果关系图简称为因果关系图。作者隐含了一个假设，即每个成分取决于其过去。【默认对任何一个 $X_a (t)$ ，过去的 $X_a$ 对其有影响】这可以用有向自环（self-loops）表示。由于这些自环的插入不改变图的分离特性（separation properties），本文省略了自环。

​	TSC图 $G_{TS}$ 和因果关系图 $G_C$ 的关系很简单：

<span id="prp2.5">**命题2.5**</span>：【聚合Aggregation】对于平稳序列 $X$ 和TSC图 $G_{TS}$ 、因果关系图 $G_C$ ，有
$$
\begin{array}{l}{(i)\qquad a \rightarrow b \notin E_{\mathrm{C}} \Leftrightarrow(a, t-u) \longrightarrow(b, t) \notin E_{\mathrm{TS}} \quad \forall u>0 \quad \forall t \in \mathbb{Z}} \\(ii)\qquad {a-b \notin E_{\mathrm{C}} \Leftrightarrow(a, t)-(b, t) \notin E_{\mathrm{TS}} \quad \forall t \in \mathbb{Z}}\end{array}
$$
**证明**：$(i)$ 根据Lauritzen（1996）的交叉性质（intersection property）之C5；$(ii)$ 显然

> 通俗解释版本：$(i)$ 的意思是， $a$ 在因果关系图中不是 $b$ 的原因等价于对任意 $t$ ，在TSC图中 $a$ 的任意过去都不是 $b$ 的原因；$(ii)$ 的意思是， $a$ 在因果关系图中和 $b$ 同期不相关等价于对任意 $t$ ，在TSC图中 $a$ 都和 $b$ 不相关

**例子2.2（延申）**：对于命题2.5和下式
$$
(a, t-u) \longrightarrow(b, t) \in E_{\mathrm{TS}} \Leftrightarrow u \in\{1, \ldots, p\} \text { and } A_{b a}(u) \neq 0
$$

可得：$a \longrightarrow b \notin E_C $ 当且仅当 $A_{ba}(1)=A_{ba}(2)=0$。得到 $X=\{X(t)\}$ 的因果图见下图

![fig2](https://raw.githubusercontent.com/rucwyf/Graphical-Model/master/pictures/Causality%20and%20graphical%20models%20in%20time%20series%20analysis/CAGMITSA2.png)

从这张图中可以看到，相对于整个过程来说，$X_1$ 与 $X_4$ 无关。更直观地说，从 1 到 4 的有向路径表面，$X_1$ 间接对 $X_4$ 造成影响。由于所有从 1 到 4 的路径都交于 3 ，这个间接原因看起来是由 $X_3$ 作为中介的。下个部分将会看到这种因果关系确实可以从图中正是推导出来。