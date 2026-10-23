# 数学公式与论文对应

本说明保留 v1.0.1 的数学定义，并说明 v1.1.0 本地候选的计算与分析方法。
下节的 `A/R` 推导继续作为 `method="quotient"` 的参考；默认内核已经改为后文的
`method="direct"`。两种方法始终保持 `H0=-1`、完整 `g[n]` 和原始 `p[n]` 的定义不变。
运行步骤与关键变量查询见 [README](../README.md)，实际验证见 [验证记录](validation.md)。

## 既有 quotient 推导与固定约定

程序首先验证原点是平衡点、线性矩阵为 `[[0,1],[-1,0]]`，并拒绝非多项式输入。
系数被视为与时间无关的常数；对真实周期的解释要求实系统。

自动生成

$$
\dot r=\sum_{j=2}^{d}G_j(\theta)r^j,\qquad
\dot\theta=H_0+\sum_{j=1}^{d-1}H_j(\theta)r^j,\qquad H_0=-1.
$$

倒数和径向展开使用

$$
A_0=\frac1{H_0},\qquad
A_n=-\frac1{H_0}\sum_{j=1}^{\min(n,d-1)}H_jA_{n-j},
\qquad R_j=\sum_{m=2}^{\min(j,d)}G_mA_{j-m}.
$$

设 `r(theta,rho)=rho+sum(u[n](theta)*rho^n,n>=2)`，逐阶提取

$$
u_n'(\theta)=[\rho^n]\sum_{j=2}^{n}R_j(\theta)
\left(\rho+\sum_{k=2}^{n-1}u_k(\theta)\rho^k\right)^j,
\qquad u_n(0)=0.
$$

程序使用截断多项式乘法，只保留需要的幂次；无需高阶求导再除以阶乘。
默认保留的中间量阶数为

$$
N_u=\max(1,N_g,N_p),\qquad N_A=\max(0,N_p,N_u-2).
$$

按用户当前约定，焦点量为

$$
g_n=\frac{u_n(2\pi)}{2\pi},\qquad n\ge2.
$$

周期部分先将径向解代入倒数展开：

$$
B_n(\theta)=-[\rho^n]\sum_{j=0}^{N_p}A_j(\theta)
r(\theta,\rho)^j,\qquad
p_n=\int_0^{2\pi}B_n(\theta)\,d\theta.
$$

因此 `B[0]=1`、`p[0]=2*Pi`，且
`P_series=2*Pi+sum(p[n]*rho^n,n=1..Np)`。
例如

$$
B_1=-A_1,\qquad B_2=-(A_1u_2+A_2),\qquad
B_3=-(A_1u_3+2A_2u_2+A_3).
$$

实际正时间对应角度递减。程序保持用户的 `0..2*Pi` 焦点量方向，只在周期积分前
加入负号。**已知系统是中心时**，该积分给出同一闭轨的正周期；对一般焦点，
输出的是这一约定下的形式时间系数。有限阶焦点量消失不足以证明中心。

所有 `g[n]` 保留完整返回映射系数，未使用低阶消失条件作约化。例如焦点示例
中的 `g[5]=a*b+3*Pi*a^2`。`Compute` 不依据有限阶消失自动认定中心、等时性或
弱中心阶数，也不自动删除奇数周期系数。后文的 `CenterCertificate` 只核验明确支持的
中心充分条件；它不构成一般中心判定算法。

## 默认 direct 内核：在原多项式方程上递推

记

$$
G(\theta,s)=\sum_{j=2}^{d}G_j(\theta)s^j,\qquad
H(\theta,s)=H_0+\sum_{j=1}^{d-1}H_j(\theta)s^j.
$$

径向方程等价于

$$
H(\theta,r)r_\theta=G(\theta,r).
$$

已经计算到第 $n-1$ 阶时，令

$$
s_{<n}=\rho+\sum_{k=2}^{n-1}u_k(\theta)\rho^k.
$$

因为 $G$ 从 $s^2$ 开始、$H-H_0$ 从 $s$ 开始，而 $u_1'=0$，未知的 $u_n$ 不会
进入下面的右端：

$$
u_n'(\theta)=\frac1{H_0}[\rho^n]\left\{
G(\theta,s_{<n})-
\bigl(H(\theta,s_{<n})-H_0\bigr)\partial_\theta s_{<n}
\right\},\qquad
u_n(\theta)=\int_0^\theta u_n'(\varphi)\,d\varphi.
$$

因此可直接用有限次数的 $G,H$ 递推，无须在径向解的复合中使用越来越多的 $R_j$。
为兼容原有查询，默认结果仍生成对应**原始向量场**的 `A`、`R` 数组；这部分成本并未
凭空消失。`quotient` 内核保留作交叉核对，两者计算同一组数学量。

周期部分先复合原始角速度：

$$
h(\theta,\rho)=H(\theta,r(\theta,\rho))
=\sum_{n\ge0}h_n(\theta)\rho^n,\qquad h_0=H_0.
$$

带方向负号的时间密度满足 $hB=-1$，从而

$$
B_0=-\frac1{H_0}=1,\qquad
B_n=-\frac1{H_0}\sum_{k=1}^{n}h_kB_{n-k},\qquad
p_n=\int_0^{2\pi}B_n(\theta)\,d\theta.
$$

在本约定 $H_0=-1$ 下，$B_n=\sum_{k=1}^n h_kB_{n-k}$。
`composed_H` 保存 direct 内核的 $h_n$，不要把它与原角速度系数 `H[n]` 混淆。

若第一个非零角速度修正为 $H_qr^q$，则 $p_N$ 至多依赖到 $u_{N-q+1}$。
`compact=true` 利用这一点缩短实际径向数组；当 $N<q$ 或角速度没有非线性修正时，
只为周期保留 $u_1$ 即可。默认 `compact=false` 继续保留旧的数组范围，实际范围始终
以 `u_order` 为准。

## 公共时间因子与次数结构

若原系统具有公共多项式因子

$$
X=q\widehat X,\qquad Y=q\widehat Y,\qquad q(0,0)=1,
$$

则在原点附近 $q$ 非零，且

$$
G=q\widehat G,\qquad H=q\widehat H,\qquad
\frac{dr}{d\theta}=\frac{\widehat G}{\widehat H}.
$$

轨道因子可在计算 $u_n,g_n$ 时约去，但时间密度仍必须使用
$-1/(q\widehat H)$。因此改变 $q$ 可以改变周期而不改变轨道返回映射。
`orbit_factor`、`orbit_G`、`orbit_H` 记录内部轨道场；公共 `G,H,A,R` 仍保持原始定义。
后面的中心示例中 $q=1+bx$，故精确径向解不含 $b$，但周期系数仍含 $b$。

程序的次数分析给出的是**可能非零阶数的上界**。设径向移位集合和角速度移位集合为

$$
S_G=\{j-1:G_j\ne0,\ j\ge2\},\qquad
S_H=\{j:H_j\ne0,\ j\ge1\},
$$

令 $M$ 为 $S_G\cup S_H$ 生成的含零加法半群。递推可出现的阶必须满足

$$
u_n\ne0\ \Longrightarrow\ n-1\in S_G+M,\qquad
B_n\ne0\ \Longrightarrow\ n\in S_H+M\quad(n\ge1).
$$

不满足必要条件的阶可以直接跳过；满足条件仍可能因系数相消或积分为零而消失。
程序不会通过浮点试算猜测结构性零。

例如仅有齐次 $d$ 次非线性时，允许的阶满足

$$
u_n,g_n:\ n\equiv1\pmod{d-1},\qquad
p_n:\ n\equiv0\pmod{d-1}.
$$

纯三次系统因而具有偶数阶 $u_n,g_n$ 为零、奇数阶 $p_n$ 为零的结构。
二次与三次混合时通常没有这一奇偶限制，不能统一删除奇数周期系数。

### 缩放律与参数权

对独立坐标缩放 $x=c\xi,y=c\eta$，相应向量场为
$F_c(\xi,\eta)=c^{-1}F(c\xi,c\eta)$，其中 $c\ne0$。
原系统 $d$ 次项的系数乘以 $c^{d-1}$，而

$$
r_c(\theta,\rho)=c^{-1}r(\theta,c\rho),\qquad
g_n^{(c)}=c^{n-1}g_n,\qquad
p_n^{(c)}=c^np_n.
$$

这提供独立的计算恒等式，也是 `ScalingCheck` 的依据。若每个原始多项式系数分别
作为独立参数，则它的权为 $d-1$，$g_n$ 的加权次数为 $n-1$，$p_n$ 的加权次数为 $n$。
同一个用户参数若同时进入不同次数的系数，不能给它强行指定一个统一权；程序按原始
各次齐次分量缩放，而不做这种错误假设。符号缩放的结果在 $c\ne0$ 的域上解释。

## 半周返回映射、完整焦点量与自然振幅

将半径视为有符号的形式变量，定义

$$
S(\rho)=r(\pi,\rho),\qquad
J(\rho)=-S(\rho),\qquad
\Pi(\rho)=r(2\pi,\rho).
$$

极坐标对 $(\theta,r)\mapsto(\theta+\pi,-r)$ 的关系给出

$$
\Pi=J\circ J.
$$

这一恒等式不要求中心成立，`HalfReturn` 返回相应的截断复合残差。
设 $S=\rho+s_2\rho^2+s_3\rho^3+\cdots$，直接复合可得

$$
[\rho^2]\Pi=0,\qquad
[\rho^3]\Pi=2(s_3-s_2^2),\qquad
[\rho^4]\Pi=s_2(s_3-s_2^2).
$$

因此对于程序的完整焦点量，

$$
g_2=0,\qquad g_4=\frac{s_2}{2}g_3.
$$

这说明某些偶数阶量依赖于低阶量，却不表示一般系统的所有偶数阶完整 $g_n$ 恒为零。
特别是 $g_3=0$ 时可推出 $g_4=0$；程序仍保留原始完整系数供查询。

### 已证中心中的偶次周期常数

若系统已经独立证明为中心，则 $\Pi=\mathrm{id}$，故 $J^2=\mathrm{id}$。
定义接近恒等的振幅变换

$$
\alpha(\rho)=\frac{\rho-J(\rho)}2,
\qquad \alpha'(0)=1.
$$

于是 $\alpha(J(\rho))=-\alpha(\rho)$。同一闭轨的周期满足 $P\circ J=P$，所以

$$
\widetilde P(\alpha)=P(\rho(\alpha))
=2\pi+\sum_{k\ge1}\tau_k\alpha^{2k}.
$$

`CenterPeriod` 同时返回振幅变换、逆变换、完整 `P_alpha_series`、`tau[k]` 和
`odd_residual`。原来的 `p[n]` 不变。默认要求支持的独立中心证书；`formal=true`
只允许做尚未证明中心的代数变换，结果会明确标成 `formal_only`，不保证奇次项消失。
此外，compact 结果可能缺少构造振幅变换所需的径向阶数，此时需要按提示补算。

对下节的中心示例，

$$
S(\rho)=\frac{\rho}{1+2a\rho},\qquad
\alpha=\frac{\rho(1+a\rho)}{1+2a\rho},\qquad
\rho=\alpha+a\alpha^2-a^3\alpha^4+O(\alpha^6).
$$

因此

$$
\tau_1=p_2=\pi b(b-a),\qquad
\tau_2=p_4-5a^2p_2
=\frac{\pi}{4}b(b-a)(3b^2-6ab-5a^2).
$$

原始 $p_3=-2ap_2$ 与这套偶次振幅展开相容，并非程序误差。
例如 $a=1,b=2$ 时，$\tau_1=2\pi$、$\tau_2=-5\pi/2$，与原始半径中的
$p_4=15\pi/2$ 是不同坐标下的系数。

临界周期的截断候选由

$$
\frac{\widetilde P'(\alpha)}{2\alpha}
=\sum_{k=1}^{K}k\tau_kz^{k-1}=0,\qquad z=\alpha^2>0
$$

给出。`CriticalPeriods` 返回候选方程、可求的根以及参数 Jacobian 的秩信息。
一般参数公式只适用于相应的非退化分支；首项系数或分母为零时必须另行分析。
截断正根和 Jacobian 的一般秩都不能单独证明真实临界周期的个数，还需要余项控制、
中心参数族及所用分岔定理的条件。

## 正物理时间的返回映射

保持当前 $g_n$ 定义时，$\Pi$ 沿增加的角度绕行；正物理时间的局部返回映射为

$$
\Pi_{\mathrm{forward}}=\Pi^{-1}.
$$

`ForwardReturn` 通过形式级数求逆，并返回复合残差。高阶系数不能逐项反号，例如

$$
\Pi(\rho)=\rho+c_3\rho^3+c_5\rho^5+\cdots
\quad\Longrightarrow\quad
\Pi^{-1}(\rho)=\rho-c_3\rho^3+(3c_3^2-c_5)\rho^5+\cdots.
$$

只有首个非零位移项的系数可以直接反号。对实系统，正时间返回图的首个可判定非零项
为负表示局部吸引、为正表示局部排斥；其前所有低阶项必须已确认为零。
参数符号未确定或只算到全零截断时，程序保留 `unknown`，不据此认定中心或稳定性。

## 独立笛卡尔同调计算

`Lyapunov(X,Y,x,y,N)` 不依赖极坐标积分。这里 $N$ 是齐次总次数上限，不是 $g$ 的下标。
写成 $X=y+\sum_{j\ge2}X_j$、$Y=-x+\sum_{j\ge2}Y_j$，构造

$$
V=x^2+y^2+\sum_{m=3}^{N}V_m,
\qquad
XV_x+YV_y=\sum_{k\ge1}L_k(x^2+y^2)^{k+1}+O(\|(x,y)\|^{N+1}).
$$

定义线性同调算子 $\mathcal L_0=y\partial_x-x\partial_y$，则每阶只需求解

$$
\mathcal L_0V_m+W_m=
\begin{cases}
0,&m\text{ 为奇数},\\
L_{m/2-1}(x^2+y^2)^{m/2},&m\text{ 为偶数},
\end{cases}
$$

其中 $W_m$ 仅由已知低阶项组成。偶数阶的径向核为
$(x^2+y^2)^{m/2}$，程序固定 $V_m$ 的 $x^m$ 系数为零以消除这一自由度。
因而这是一套明确归一化的独立 Lyapunov 量；`residual=0` 表示所算总次数以内的
多项式恒等式成立，不表示完整无限级数恒等式已被证明。

在 $L_1=\cdots=L_{k-1}=0$、且本程序归一化和方向约定成立的条件下，首个可能非零量满足

$$
g_{2k+1}=-\frac{L_k}{2}.
$$

原因是可取形式半径 $R=\sqrt V=r+O(r^2)$，其首个漂移项在物理时间中为
$\dot R=(L_k/2)R^{2k+1}+\cdots$，而 $\dot\theta=-1+\cdots$。
截面上的 $R=\rho+O(\rho^2)$ 是接近恒等的变换，不改变首个非零返回位移的系数；
沿 $0\to2\pi$ 积分后再除以 $2\pi$，得到上述负号与因子。
低阶 $L$ 不消失时，高阶完整 $g_n$ 含额外混合项，不能把所有 `L[k]` 逐项乘以
$-1/2$ 充当原始 `g[2*k+1]`。输出有限多个零 `L` 同样不足以证明中心。

## 参数局部展开、条件域与中心证据

### 参数 jet

`ParameterJet` 沿指定方向代入

$$
\lambda=\lambda_0+\varepsilon v
$$

并在系数环 $\mathbb K[\varepsilon]/(\varepsilon^{K+1})$ 中推进递推。
截断贯穿乘积、化简和积分，不是先求完整参数通式后再删高次项。每个方向参数都必须
给出基点值；状态变量不能作为方向参数，基点或方向也不能依赖状态变量。
基点分母为零等导致的奇异展开会被拒绝。

返回的 `g_coefficients[n,k]`、`p_coefficients[n,k]` 是 $\varepsilon^k$ 的 **Taylor 系数**；
若需要第 $k$ 阶导数，须再乘 $k!$。私有展开变量通过 `ans["epsilon"]` 取得，不与用户
模型中已有的同名参数混淆。局部方向展开不等于完整参数通式，方向若离开中心参数族，
相应 `p_n` 仍只能解释为形式时间系数。

### 多项式理想与非零条件

`ZeroStatus` 与 `AnalyzeConditions` 分别保留原式、约化式、方程关系、非零前提和分母域。
对有理式先在已记录的分母非零域上约化分子，使用 Gröbner 基求正规余式。
余式为零给出理想成员的充分依据；余式非零一般只返回 `unknown`，因为它仍可能在
某些允许参数点上取零。明确的非零常数或用户给定的非零条件才可提供相应非零依据。
此处使用 Maple 的 [`Groebner:-NormalForm`](https://www.maplesoft.com/support/help/maple/view.aspx?path=Groebner/NormalForm)
作为代数运算基础。

例如 $a(b-a)=0$ 同时包含 $a=0$ 与 $b=a$ 两个分支，不能擅自约掉 $a$。
当前分析不声称完成根理想成员判定、非零条件饱和或全部实分支求解；`unknown` 也不代表
数学上不可判定。实参数可行性与复代数理想关系必须分别解释。
`AnalyzeConditions` 将原始输入的分母非零条件一并传入，防止某个已相消为零的系数掩盖
原系统的不可定义参数点。分母记录仅覆盖当前表达式仍可见的部分，输入形成前已经
被 Maple 约掉的因子无法恢复。

### 独立中心充分条件

`CenterCertificate` 重新核验原点平衡与固定线性部分，在明确列出的实参数及分母域前提下，
支持以下充分条件：

- $X_x+Y_y=0$：局部 Hamilton 函数的二次部分为 $(x^2+y^2)/2$，给出中心。
- $X(x,-y)=-X(x,y)$、$Y(x,-y)=Y(x,y)$：实解析系统关于 $y\mapsto-y$ 可逆，配合固定非退化椭圆线性部分给出中心。
- 给定实多项式第一积分 $I$，核验 $I_xX+I_yY=0$、$\nabla I(0)=0$，且 Hessian 正定或负定。

未知参数按实参数前提明确记录，未能核验实系数、浮点输入、未核验的一般解析第一积分、
未能判定定号等情况返回 `unknown`。`proved` 表示所列前提下已通过对应充分条件，
不会抹去这些前提，也不由有限阶焦点量消失推出。

`InputInfo` 和 `RunInfo` 保存原输入、计算时环境、阶数及方向约定，便于比较结果。
显式给出 `source_file` 时可以计算文件 SHA1 指纹，用于区分文件；该指纹既不是安全签名，
也不能自动证明指定文件就是已经加载的源码。

## 示例结果

解析中心示例为

$$
\dot x=(1+bx)(y+axy),\qquad
\dot y=(1+bx)(-x+ay^2).
$$

它给出

$$
G_2=a\sin\theta,\quad G_3=ab\cos\theta\sin\theta,\quad
H_0=-1,\quad H_1=-b\cos\theta,\quad H_2=0.
$$

精确径向解是

$$
r(\theta,\rho)=\frac{\rho}{1+a\rho(1-\cos\theta)}.
$$

该解对小初始半径为正且以 `2*Pi` 为周期，所以所有焦点量为零。计算得到

$$
p_1=0,\quad p_2=\pi b(b-a),\quad
p_3=-2\pi ab(b-a),\quad
p_4=\frac{3\pi}{4}b(b-a)(b^2-2ab+5a^2).
$$

取 `a=1,b=2`：

$$
P(\rho)=2\pi+2\pi\rho^2-4\pi\rho^3+
\frac{15\pi}{2}\rho^4+O(\rho^5).
$$

另一个径向三次焦点示例验证 `g[3]=-a` 和
`g[5]=a*b+3*Pi*a^2`；其形式 `p[4]` 包含漂移项 `4*Pi^2*a*b`。

## 与论文的对应及验证

依据用户提供的 Zhang, Hou, Zeng，*Weak centers and bifurcation of critical
periods in reversible cubic systems*，*Computers & Mathematics with
Applications* **40** (2000), 771–782，第 2 节式 (6)–(15)，以及第 781–782 页附录。

原文系统的线性部分为 `[[0,-1],[1,0]]`，给出 `H0=1`。本程序使用包含 `H0` 的统一
公式，适配用户的反向线性旋转，并推广到任意多项式次数。焦点量的编号及归一化
采用用户明确给出的定义；论文的重点是中心的周期系数。

核对论文时，将原文式 (1) 的**整个向量场反号**，使输入符合用户标准。该操作
只用于这个验证示例，程序不会自动反号输入。闭轨及其正周期保持不变。
论文第 774 页公布的完整公式为

$$
p_2=\frac{\pi}{12}\left(
4a_1^2+10a_1a_2+10a_2^2-5a_1b_1-a_2b_1+b_1^2
+3a_3+9a_4-9b_2-3b_3\right).
$$

实际计算与该公式完全一致。

### 新增方法的来源说明

本轮 direct 递推、公共时间因子分离、次数支撑与缩放恒等式、半周复合和自然振幅公式，
均从本文固定的径向方程与形式级数关系自行推导；没有复制第三方代码。
Gröbner 约化调用 Maple 内置库，具体正规余式含义参见上述官方文档。

独立笛卡尔计算参考了直接在直角坐标中计算 Poincaré–Lyapunov 中心条件这一研究方向：
[Boularas（2009），*Théorème de Poincaré–Lyapunov et conditions nécessaires d'existence d'un centre*](https://comptes-rendus.academie-sciences.fr/mathematique/articles/10.1016/j.crma.2009.07.009/)。
本程序的同调矩阵、$V_m$ 归一化及与当前 $g_n$ 的首非零项关系在本文单独写明，
不声称与该文的全部算法或编号相同，也没有改编其代码。
