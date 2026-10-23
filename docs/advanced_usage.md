# v1.1.0 候选高级用法

本文对应 **v1.1.0 开发候选，尚未发布**。现有公开发行版及其单文件 `focus_period.mpl` 仍为 v1.0.1，不能运行本文新增的接口或关键字。所有新功能仍在同一个 `src/focus_period.mpl` 中，只需加载该文件；没有新增运行依赖。

路径、首次加载、旧变量及全部下标查询见 [README](../README.md)。实际验证范围见[验证记录](validation.md)，公式与方向约定见[数学说明](mathematics.md)。发布新版本须另外取得用户确认。

## 1. 加载与接口速查

以下代码从项目根目录运行。Windows 字符串使用 `/`，或将每个反斜杠写成 `\\`；不要原样粘贴文件管理器的单反斜杠路径。在 Maple 的一维可执行输入区运行即可，无需改成二维公式。

```maple
restart:
currentdir("D:/Maple/Maple_Focus_Period"):
read "src/focus_period.mpl":
FocusPeriod:-Version();  # 本文候选源码应显示 "1.1.0"。

X := (1+b*x)*(y+a*x*y):
Y := (1+b*x)*(-x+a*y^2):
ans := FocusPeriod:-Compute(X,Y,x,y,5,4):
```

目录或 `read` 报错时，先修正第一个路径错误，再调用模块。更新磁盘文件后，需要重新 `restart`、加载与定义系统；保存文件本身不会更新已加载模块。

| 接口 | 主要用途与返回内容 |
|---|---|
| `Compute(X,Y,x,y,Ng,Np,...)` | 原六参数调用保留；生成原始 g/p 和全部旧中间量，关键字控制计算内核 |
| `PolarData(X,Y,x,y)` | 仅提取 G/H，不递推 u/g/p；也接受 `jet=[epsilon,K]` |
| `InputInfo(X,Y,x,y)` | 原始输入、参数、精确性、可见分母、环境与数学约定 |
| `ZeroStatus(f,...)` | 在指定代数条件下给出 `zero/nonzero/unknown` 及依据 |
| `AnalyzeConditions(ans,relations,parameters,...)` | 分别分析已算 g/p，保留原始结果不变 |
| `CenterCertificate(ans,...)` | 核验少数受支持的实中心充分条件；不由有限 g 消失推断中心 |
| `RunInfo(ans,...)` | 记录环境、阶数及可选显式源文件指纹 |
| `ParameterJet(X,Y,x,y,Ng,Np,base,direction,K)` | 沿指定参数方向，直接计算到 K 阶局部参数展开 |
| `HalfReturn(ans)` | 半周图、带符号半周图及全周复合残差 |
| `ForwardReturn(ans)` | 正物理时间的形式返回图及有充分符号依据时的稳定性判断 |
| `CenterPeriod(ans,...)` | 自然振幅下的偶次周期常数；默认要求独立中心证据 |
| `CriticalPeriods(period_data,parameters)` | 截断临界周期候选与周期常数 Jacobian 的 generic 秩 |
| `Lyapunov(X,Y,x,y,N)` | 独立笛卡尔坐标内核；N 是最高总次数，所得 L 与原始 g 区分 |
| `ScalingCheck(ans,scale)` | 独立重算缩放后的系统，检查系数的缩放权 |

`Show(ans)` 仍用于简洁展示基础计算结果。新增分析接口均显式调用，避免每次 `Compute` 自动重复进行所有研究分析。

完整可执行演示见 [examples/04_research.mpl](../examples/04_research.mpl)。全部本地Maple回归可用 `tools/verify.ps1` 运行；冷启动性能比较用 `tools/benchmark.ps1`，具体命令和测量边界见[验证记录](validation.md)。

## 2. 计算内核、次数支撑与共同因子

### 两条内核如何选择

```maple
direct_result := FocusPeriod:-Compute(X,Y,x,y,5,4):
quotient_result := FocusPeriod:-Compute(X,Y,x,y,5,4,
    method="quotient"):

simplify(direct_result["g"][5]-quotient_result["g"][5]);
simplify(direct_result["p"][4]-quotient_result["p"][4]);
# 对该示例均为 0。完整核验应比较各目标阶及方程残差。
```

默认 `method="direct"` 从

\[
H(\theta,r)\,r_\theta=G(\theta,r)
\]

直接递推径向解。时间部分先将径向解代入 H，再解 `H*B=-1`，因此仍有 `B[0]=1`、`p[0]=2*Pi`。

`method="quotient"` 保留先求 A/R、再复合径向解和时间密度的原路径。两者的数学目标及原始 g/p 定义相同。不同模型的表达式增长方式不同，不能保证默认路径对所有输入都有相同提速。

### 默认值与保留范围

```maple
ans := FocusPeriod:-Compute(X,Y,x,y,5,4,
    method="direct", optimize=true, factor_orbits=true,
    compact=false, jet=[]):

ans["method"]; ans["optimize"]; ans["compact"]; ans["jet"];
ans["u_order"]; ans["A_order"];
```

默认 `compact=false` 保留旧范围：

\[
N_u=\max(1,N_g,N_p),\qquad N_A=\max(0,N_p,N_u-2).
\]

全部旧结果字段与查询方式保留。`optimize=true` 只跳过能够根据次数结构证明为零的阶，数组内仍填入零；`optimize=false` 使用逐阶稠密递推。

```maple
dense_result := FocusPeriod:-Compute(X,Y,x,y,5,4,optimize=false):
ans["structure"]["radial_shifts"];
ans["structure"]["angular_shifts"];
ans["orbit_structure"]["u_possible_orders"];
ans["structure"]["period_possible_orders"];
```

`possible_orders` 是次数允许出现的阶数上界。列出的阶仍可能因积分或系数抵消而为零；未列出的阶才是由该结构分析排除的阶。

只有显式设置 `compact=true` 才缩减周期计算需要的 u 阶数，仍会满足请求的 Ng。查询实际范围时使用 `u_order`，不要继续假定等于 `max(1,Ng,Np)`。

```maple
circular := FocusPeriod:-Compute(y*(1+x^2+y^2),
    -x*(1+x^2+y^2),x,y,0,4,compact=true):
circular["u_order"];            # 3，默认非 compact 时为 4。
circular["period_u_order"];     # 本次周期计算所需的径向阶数。
circular["requested_u_order"];  # 旧默认范围对应的阶数。
```

自然振幅转换可能需要比计算 p 更高的 u 阶数。准备调用 `CenterPeriod(ans)` 时，保留默认 `compact=false`，或让 Ng 至少达到 Np；缺少阶数时会明确报错。

### 轨道因子与原始场分开保存

共同因子 q 归一化为 `q(0,0)=1`。它在原点附近非零时可以从径向轨道方程约去，但仍须保留在时间计算中。

```maple
ans["orbit_factor"];             # 本例为 1+b*x。
ans["orbit_degree"];             # 本例约去因子后为 2。
ans["orbit_G"]; ans["orbit_H"];  # 简化轨道的极坐标数据。
ans["G"]; ans["H"];             # 原始系统的数据。
ans["A"]; ans["R"];             # 始终按原始场定义。

without_factor := FocusPeriod:-Compute(X,Y,x,y,5,4,
    factor_orbits=false):
simplify(ans["p"][4]-without_factor["p"][4]);  # 0。
```

本例径向解与 b 无关，周期仍依赖 b。不能通过删除共同因子来代替真实周期计算。

## 3. 参数条件、分母与三态零判定

```maple
input_info := FocusPeriod:-InputInfo(X,Y,x,y):
input_info["parameters"];
input_info["exactness"];             # "exact" 或 "floating"。
input_info["denominator_nonzero"];
input_info["real_parameter_assumption"];
# Compute 的记录也可通过 ans["input"] 访问。
```

`1/10` 是精确数，`0.1` 会产生浮点输入。浮点近零不作为精确代数零证明。分母记录只覆盖表达式当前仍可见的分母；输入形成前已经被 Maple 约掉的因子无法恢复。

```maple
unknown_case := FocusPeriod:-ZeroStatus(a*(b-a),parameters=[a,b]):
zero_case := FocusPeriod:-ZeroStatus(a*(b-a),
    relations=[a=0],parameters=[a,b]):
nonzero_case := FocusPeriod:-ZeroStatus(b,
    parameters=[a,b],nonzero=[b]):
unknown_case["status"];  # "unknown"。
zero_case["status"];     # "zero"，有条件 a=0。
nonzero_case["status"];  # "nonzero"，使用显式前提 b<>0。

conditions := FocusPeriod:-AnalyzeConditions(ans,[a=0],[a,b]):
conditions["g"][3]["status"];
conditions["p"][2]["reduced"];
conditions["p"][2]["evidence"];
conditions["p"][2]["denominator_nonzero"];
```

`relations` 是等于零的多项式关系或等式列表；`parameters` 必须列出关系中所有符号参数。`nonzero=[b]` 表示 b 必须非零。

条件约化使用多项式理想。余式为零可提供相应理想条件下的恒等式依据；余式非零通常仍可在特殊参数点消失，所以状态保持 `unknown`。本功能不自动计算根理想、不穷尽实分支、不做非零条件饱和，也不证明参数域非空。矛盾条件记录在 `domain_status` 中。

例如 `a*(b-a)=0` 不能直接除以 a；`a=0` 与 `b=a` 应分别分析。`AnalyzeConditions` 返回独立表，不改写 `ans["g"]` 或 `ans["p"]`。

## 4. 中心充分条件与自然振幅周期

### 中心核验

```maple
certificate := FocusPeriod:-CenterCertificate(ans):
certificate["status"];       # 本例为 "proved"。
certificate["method"];       # 本例有关于 y 反射的可逆性。
certificate["assumptions"];  # 参数实性等适用前提。
certificate["evidence"];
```

当前覆盖以下充分条件：实解析场的零散度与固定非退化椭圆线性部分；关于 `(x,y)->(x,-y)` 的可逆性；或给定多项式第一积分满足恒等式且在原点有非退化极值。

```maple
integrable := FocusPeriod:-Compute((1+x+y)*y,-(1+x+y)*x,x,y,4,4):
certificate_I := FocusPeriod:-CenterCertificate(integrable,
    first_integral=x^2+y^2):
certificate_I["status"];
certificate_I["first_integral_residual"];
```

未知参数按记录中的实参数前提处理；参数分式还须满足输入记录中的分母非零条件，不能将奇异参数点纳入结论。复系数、浮点数据或无法核实的实性不会获得该中心证明。检查失败返回 `unknown`，不能解释为已证明存在焦点。若更简单的充分条件已成立，过程可能直接采用该条件，不再检查用户提供的第一积分。

### 半周图与原始 g 的依赖

```maple
half := FocusPeriod:-HalfReturn(ans):
half["S_series"];                # S(rho)=r(Pi,rho)。
half["J_series"];                # J=-S。
half["full_return_series"];
half["composition_residual"];    # J(J(rho))-r(2*Pi,rho) 的截断残差。
half["dependencies"]["g4_minus_s2_g3_over_2"];
```

带符号半周图满足全周复合关系。若 `s2=[rho^2]S`，则 `g[2]=0`、`g[4]=s2*g[3]/2`。一般混合次数系统的完整偶数阶 g 不必全部为零。截断复合一致也不构成中心证明。

### 自然振幅与独立偶次周期常数

```maple
period_data := FocusPeriod:-CenterPeriod(ans):
al := period_data["alpha"]:
period_data["alpha_series"];      # alpha(rho)=(rho-J(rho))/2。
period_data["rho_of_alpha"];      # 反变换，使用返回的私有 alpha。
period_data["P_alpha_series"];
period_data["tau"][1];           # alpha^2 系数。
period_data["tau"][2];           # alpha^4 系数。
period_data["odd_residual"];     # 本例为 0。
period_data["center_status"];
```

在已证明的中心上，新振幅将带符号半周图化为 `alpha->-alpha`，周期成为

\[
\widetilde P(\alpha)=2\pi+\tau_1\alpha^2+\tau_2\alpha^4+\cdots.
\]

原始 rho 坐标下的奇数周期项可以非零；本例有 `p[3]=-2*a*p[2]`。新 tau 单独保存，不删除或改写原始 p。

可用 `CenterPeriod(ans,N)` 指定不高于已算 Np 的转换阶数。默认要求 `CenterCertificate` 成功，也可传入 `first_integral=...`。`formal=true` 允许在没有证明时进行代数转换，但返回 `center_status="formal_only"`；它不证明中心，也不保证奇次项消失。

### 截断临界周期候选与参数秩

```maple
critical := FocusPeriod:-CriticalPeriods(period_data,[a,b]):
critical["critical_polynomial"];
critical["candidate_domain"];      # z>0，z=alpha^2。
critical["candidates"];
critical["generic_branch_conditions"];
critical["jacobian"];
critical["generic_rank"];
```

候选方程为

\[
\tau_1+2\tau_2z+3\tau_3z^2+\cdots=0,\qquad z=\alpha^2>0.
\]

当前对不高于四次的非零截断候选多项式调用符号求解；更高次数保留方程。候选的正性依赖参数时标为未确定；首项退化、分母为零及其他特殊参数分支要另行处理。

只有 `odd_residual=0` 时才能用 `z=alpha^2` 表示完整截断导数；含奇次项的形式转换会被 `CriticalPeriods` 拒绝，避免只分析偶次投影而丢项。

Jacobian 的 `generic_rank` 是一般参数位置的代数秩，特殊分支可能降低。参数变化还必须保留已确认的中心族。截断正根及满秩都不能单独证明真实临界周期的存在、个数或局部分岔上界；真实结论仍需余项控制和相应定理。

## 5. 局部参数展开与方向导数

只关心某参数点附近时，无需先生成完整多参数通式。以下沿 `(a,b)=(1,2)+epsilon*(1,0)` 展开一阶：

```maple
jet_result := FocusPeriod:-ParameterJet(X,Y,x,y,5,4,
    {a=1,b=2},{a=1,b=0},1):
ep := jet_result["epsilon"]:
jet_result["p"][2];    # 2*Pi-2*Pi*ep。
jet_result["parameter_jet"]["p_coefficients"][2,0]; # 基点 p2=2*Pi。
jet_result["parameter_jet"]["p_coefficients"][2,1]; # 方向导数 -2*Pi。
jet_result["parameter_jet"]["g_coefficients"][3,1];
jet_result["parameter_jet"]["base"];
jet_result["parameter_jet"]["direction"];
```

`base`、`direction` 使用集合中的 `参数=值`；方向中每个参数都必须有基点值，未指定方向的基点参数保持固定。未列入基点的其他参数仍保留符号。

对于 K 阶展开，`p_coefficients[n,k]` 是 epsilon 的 k 次 Taylor 系数，等于 k 阶方向导数除以 `k!`；g 同理。它不是全参数通式，也不是自动求出的多参数偏导张量。

底层 `Compute(...,jet=[ep,K])` 在中间乘积、截断、化简、积分递推中持续限制 epsilon 阶数。参数分式在基点有普通 Taylor 展开时可用；有极点的基点会报错，不能将 Laurent 主部丢弃后当作正常展开。

方向扰动不一定停留在中心族内。离开中心族时，p 仍按程序约定解释为形式时间系数。需要物理周期的导数时，应先采用中心族内的参数坐标。

## 6. 正时间返回图与独立笛卡尔量

### 正时间方向

```maple
Xf := y+a*x*(x^2+y^2)+b*y*(x^2+y^2):
Yf := -x+a*y*(x^2+y^2)-b*x*(x^2+y^2):
focus := FocusPeriod:-Compute(Xf,Yf,x,y,5,4):
forward := FocusPeriod:-ForwardReturn(focus):
focus["g"][3];       # -a。
forward["g"][3];     # a。
forward["g"][5];     # 3*Pi*a^2-a*b，不能简单取 -focus["g"][5]。
forward["return_series"];
forward["inverse_residual"];
forward["stability"];
```

原始图沿 theta 增大，正时间图是它的形式逆。`ForwardReturn` 不更改极坐标符号或原始 g。只在首个非零项具有足够的精确实符号依据时判断稳定或不稳定；一般含参数时保持 `unknown`，不能把有限阶零值解释为中心。

### 笛卡尔坐标独立内核

```maple
cart := FocusPeriod:-Lyapunov(Xf,Yf,x,y,6):
cart["V"];
cart["V_components"];
eval(cart["L"]);
cart["L"][1];       # 该例为 2*a。
cart["residual"];   # 所请求总次数内，所构造恒等式的残差。
cart["normalization"];
```

这里 `N=6` 指 x/y 的最高总次数，不是 `g[6]` 或 rho 的阶数。计算构造

\[
V=x^2+y^2+V_3+\cdots+V_N,\qquad
X V_x+Y V_y=\sum_k L_k(x^2+y^2)^{k+1}
\]

至该总次数；偶数次 `V_m` 的 `x^m` 系数固定为零以选择归一化。`L[k]` 存于 table，范围为 `1..floor(N/2)-1`；`N=2` 或 `3` 时无 L 项。

L 与原始完整 g 不同，不能逐项替换。在 `L[1]..L[k-1]` 已消失的前提下，本归一化与当前角向约定给出 `g[2*k+1]=-L[k]/2`。残差零验证构造的一致性，本身不证明所有高阶阻碍都消失。

## 7. 缩放核对与可复现记录

```maple
scaling := FocusPeriod:-ScalingCheck(ans,2):
scaling["status"];
eval(scaling["g_residual"]);
eval(scaling["p_residual"]);
```

该接口独立重算坐标缩放后的系统，检查 g 的权 `n-1` 和 p 的权 `n`。缩放因子须非零；浮点输入只作为近似检查。它不会在每次 `Compute` 内自动重复计算，也不会将参数名称随意赋予权重。代数缩放恒等式与物理初始半径的解释应分开处理。

```maple
run_info := FocusPeriod:-RunInfo(ans,
    source_file="src/focus_period.mpl",revision="local-candidate"):
run_info["program_version"];
run_info["kernel_version"];
run_info["orders"];
run_info["input"];
run_info["source_digest_algorithm"]; # "sha1"。
run_info["source_digest"];
run_info["source_identity_scope"];
```

`source_file` 是显式指定的文件，可省略；省略后不生成源文件指纹。Maple 内置 `FileTools:-Hash` 在本实现使用 SHA1，仅用于辨识文件版本。它不会自动发现本次 `read` 加载了哪一个路径，指定文件的指纹也不能证明内存模块来自该文件。发行候选材料如需 SHA256，应在外部另行计算。

复现研究计算时，保留原始系统、参数关系与非零条件、精确或浮点输入、阶数、内核选项、Maple 环境和确切源文件。不同条件的结果不要仅因版本号相同而视为同一次计算。

## 8. 方法来源与验证边界

- 极坐标方向与 g/p 定义沿用本项目用户约定，参考 Zhang、Hou、Zeng（2000）的递推框架；原文与当前 H0 符号区别见[数学说明](mathematics.md)。
- 交叉相乘递推、次数支撑、共同因子分离、半周复合、自然振幅及正时间求逆根据这些方程关系推导；没有复制外部程序代码。
- 条件约化与笛卡尔线性系统使用 Maple 内置 [Groebner:-NormalForm](https://www.maplesoft.com/support/help/maple/view.aspx?path=Groebner/NormalForm) 和 [LinearAlgebra:-LinearSolve](https://www.maplesoft.com/support/help/maple/view.aspx?path=LinearAlgebra/LinearSolve)。文件辨识使用 [FileTools:-Hash](https://www.maplesoft.com/support/help/maple/view.aspx?path=FileTools/Hash)，没有引入额外库。
- 独立笛卡尔方法的方向参考 Boularas（2009），*Théorème de Poincaré–Lyapunov et conditions nécessaires d'existence d'un centre*，[出版页面](https://comptes-rendus.academie-sciences.fr/mathematique/articles/10.1016/j.crma.2009.07.009/)。本实现的具体归一化、输出与核对关系由仓库代码和测试明确给出。
- 专项与集成验证结果见[验证记录](validation.md)。有限阶 g 消失、截断残差零、形式变换、候选根和 generic 秩各有自己的证据范围，不能相互替代为完整中心或分岔证明。
