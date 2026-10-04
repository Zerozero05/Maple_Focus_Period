# 数学公式与论文对应

本说明对应 v1.0.0。运行步骤见 [README](../README.md)，实际验证见 [验证记录](validation.md)。

## 计算公式

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
计算需要的阶数为

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
中的 `g[5]=a*b+3*Pi*a^2`。程序不自动证明中心、等时性或弱中心阶数，也不自动
删除奇数周期系数。

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
