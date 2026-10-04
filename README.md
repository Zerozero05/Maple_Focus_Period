# Maple_Focus_Period

从原始标准平面多项式微分方程组出发，自动计算极坐标系数、径向递推系数、焦点量和周期系数的 Maple 程序。主程序带有中文注释，只依赖 Maple 内置功能。

当前版本：**v1.0.0**。程序已在 **Maple 2025.1** 中实际运行验证。其他 Maple 版本尚未逐一测试。

[GitHub 仓库](https://github.com/Zerozero05/Maple_Focus_Period) · [数学说明](docs/mathematics.md) · [验证记录](docs/validation.md) · [版本变更](CHANGELOG.md)

## 1. 这个程序接收什么，输出什么？

输入为微分方程右端的两个多项式表达式：

$$
\dot x=X(x,y),\qquad \dot y=Y(x,y),\qquad
D(X,Y)(0,0)=\begin{pmatrix}0&1\\-1&0\end{pmatrix}.
$$

例如，输入的是 `X := y+a*x*(x^2+y^2)`，而不是 `diff(x(t),t)=...`。原点必须是平衡点；多项式系数作为与时间无关的常数处理。程序支持任意多项式次数，计算阶数由使用者选择。

计算流程为：

```text
原始多项式 X、Y
    → 极坐标化，提取 G[i]、H[i]
    → 角速度倒数系数 A[i]
    → 径向方程系数 R[i]，径向解系数 u[i]
    → 焦点量 g[n] 和周期系数 p[n]
```

程序始终保留

$$
x=r\cos\theta,\qquad y=r\sin\theta,\qquad H_0=-1.
$$

焦点量按当前约定计算：

$$
g_n=\frac{u_n(2\pi)}{2\pi},\qquad n\ge 2.
$$

周期部分取

$$
P(\rho)=-\int_0^{2\pi}
\frac{d\theta}{H(\theta,r(\theta,\rho))}
=2\pi+\sum_{n\ge1}p_n\rho^n,
$$

因此 `P(0)=2*Pi`。**已知系统是中心时，这给出闭轨的正周期；对于一般焦点，输出的是上述约定下的形式时间系数。** 有限阶焦点量为零不足以单独证明中心。

## 2. 下载到电脑

首次使用建议下载整个项目，这样示例、主程序和工作表的相对路径都能保持一致。

1. 打开 [GitHub 仓库](https://github.com/Zerozero05/Maple_Focus_Period)。
2. 获取固定版本时，打开 [v1.0.0 发布页](https://github.com/Zerozero05/Maple_Focus_Period/releases/tag/v1.0.0)，下载 **Source code (zip)**；获取当前开发版时，点击 **Code → Download ZIP**。
3. 解压。版本包的文件夹通常名为 `Maple_Focus_Period-1.0.0`，当前开发版通常名为 `Maple_Focus_Period-main`。
4. 将解压后的项目文件夹改名为 `Maple_Focus_Period`，放到自己方便使用的位置。

以下说明统一假设项目保存在：

```text
D:/Maple/Maple_Focus_Period
```

这只是示例位置。请换成自己的实际路径。该目录内应当直接有 `src`、`examples`、`tests` 和 `docs` 四个文件夹，而不是还隔着一层 ZIP 解压目录。

## 3. `.mpl` 和 `.mw` 有什么区别？

| 文件 | 用途 | 如何使用 |
|---|---|---|
| [src/focus_period.mpl](src/focus_period.mpl) | 主程序，定义 `FocusPeriod` 模块 | 在 Maple 中用 `read` 加载；单独加载不会计算示例 |
| [examples/01_center.mpl](examples/01_center.mpl) | 已知中心的计算示例 | 在 Maple 中用 `read` 执行 |
| [examples/02_focus.mpl](examples/02_focus.mpl) | 焦点量及符号方向的计算示例 | 在 Maple 中用 `read` 执行 |
| [examples/user_system.mpl](examples/user_system.mpl) | 自己的系统的编辑入口 | 修改 `X`、`Y`、`Ng`、`Np` 后执行 |
| [examples/FocusPeriod_Worksheet.mw](examples/FocusPeriod_Worksheet.mw) | Maple 工作表入口 | 用 Maple 打开，设置路径后按工作表中的步骤运行 |

`.mpl` 是可读的 Maple 程序文本；`.mw` 是 Maple 工作表。工作表仍然读取同一个主程序，所以使用 `.mw` 时也要保留项目中的 `src` 文件夹。

**不要仅双击 `.mpl`，然后期待自动出现计算结果。** 双击可能只是打开文本编辑器；计算需要在 Maple 中执行相应的 `read` 语句。第一次使用时，可以按下一节直接运行一个已验证的示例。

## 4. 第一次运行：先算一个中心示例

打开 Maple，将下面的语句放入可执行输入区，并按顺序执行：

```maple
restart:
currentdir("D:/Maple/Maple_Focus_Period"):
read "examples/01_center.mpl";
```

`currentdir` 设置 Maple 寻找相对路径文件时使用的目录。每个新会话先设置一次项目根目录，之后就可以按本说明中的相对路径读取文件。

该示例使用

$$
X=(1+bx)(y+axy),\qquad Y=(1+bx)(-x+ay^2).
$$

预期中间量为

$$
G_2=a\sin\theta,\quad G_3=ab\cos\theta\sin\theta,\quad
H_0=-1,\quad H_1=-b\cos\theta,\quad H_2=0.
$$

其精确径向解为

$$
r(\theta,\rho)=\frac{\rho}{1+a\rho(1-\cos\theta)},
$$

因而所有焦点量都为零。低阶周期系数为

$$
p_1=0,\qquad p_2=\pi b(b-a),\qquad
p_3=-2\pi ab(b-a),
$$

$$
p_4=\frac{3\pi}{4}b(b-a)(b^2-2ab+5a^2).
$$

取 `a=1,b=2` 时，预期周期展开为

$$
P(\rho)=2\pi+2\pi\rho^2-4\pi\rho^3+
\frac{15\pi}{2}\rho^4+O(\rho^5).
$$

焦点示例可在项目根目录中执行：

```maple
read "examples/02_focus.mpl";
```

其中 `g[3]=-a`、`g[5]=a*b+3*Pi*a^2`，可以检查程序保留的角度方向和归一化约定。

### 在 `.mw` 工作表中运行

1. 保留整个项目文件夹，在 Maple 中打开 `examples/FocusPeriod_Worksheet.mw`。
2. 在第一个可执行输入块中，将 `repoRoot` 改为自己的项目根目录，例如 `"D:/Maple/Maple_Focus_Period"`。
3. 从上到下执行各输入块，先加载程序，再定义系统和计算阶数，最后查看结果。
4. 分析自己的模型时，修改工作表中的 `X`、`Y`、`Ng`、`Np`；建议另存到项目根目录下自建的 `models` 文件夹。
5. 工作表默认使用已验证的中心示例。核心算法始终来自 `src/focus_period.mpl`，不需要把算法代码粘进工作表。

## 5. 改成自己的系统

最方便的方式是复制 [examples/user_system.mpl](examples/user_system.mpl) 并改写其中的 `X`、`Y`、`Ng`、`Np`。也可以在工作表中直接运行：

```maple
restart:
currentdir("D:/Maple/Maple_Focus_Period"):
read "src/focus_period.mpl":

# 替换成你的原始系统右端。
X := (1+b*x)*(y+a*x*y):
Y := (1+b*x)*(-x+a*y^2):

# 焦点量与周期系数的最高阶，分别指定。
Ng := 5:
Np := 4:

ans := FocusPeriod:-Compute(X, Y, x, y, Ng, Np):
FocusPeriod:-Show(ans);
```

这里 `Ng` 和 `Np` 都是初始半径 `rho` 的幂次上限：

- `Ng=5`：计算并返回 `g[2]` 至 `g[5]`。
- `Ng=0`：不输出焦点量。
- `Np=4`：计算并返回 `p[1]` 至 `p[4]`；`p[0]` 始终为 `2*Pi`。
- `Np=0`：不计算非常数周期项，只保留 `P(0)=2*Pi`。
- `Ng=1` 不接受，因为 `u[1]=1` 是径向展开的初始项，不属于返回位移。

`Ng`、`Np` 必须是非负整数。只算一种量时，关闭另一种输出即可：

```maple
# 只算焦点量，到 rho^7。
only_g := FocusPeriod:-Compute(X, Y, x, y, 7, 0):
FocusPeriod:-Show(only_g);

# 只算周期系数，到 rho^6。
only_p := FocusPeriod:-Compute(X, Y, x, y, 0, 6):
FocusPeriod:-Show(only_p);
```

只算周期系数仍然需要部分径向解系数 `u[n]`；程序会自动计算所需中间量。

### 参数代入和参数假设

若只需要一组参数的结果，建议**先把参数值代入系统，再计算**。这样通常比先生成多参数高阶通式更快：

```maple
# X、Y 先按上面的方式定义为含参数 a、b 的表达式。
parameter_values := {a=1, b=2}:
X_value := eval(X, parameter_values):
Y_value := eval(Y, parameter_values):

ans_value := FocusPeriod:-Compute(X_value, Y_value, x, y, 5, 4):
FocusPeriod:-Show(ans_value);
```

使用集合将一组参数集中列出，方便检查是否遗漏。精确计算优先用 `1/3` 等有理数；`0.333333` 是近似数，可能使应当精确为零的表达式留下舍入误差。

已经算出符号结果时，也可以集中代入：

```maple
eval(ans["p"][4], {a=1, b=2});
eval(ans["P_series"], {a=1, b=2});
```

需要附加实性或不等式条件时，可以对具体的后续化简使用 `assuming`，将这一组假设写在同一条语句中：

```maple
simplify(ans["p"][4]) assuming a::real, b::real;
simplify(ans["p"][2]) assuming a>0, b>a;
```

这些假设作用于对应的化简语句，不会改写系统，也不能替代中心条件的数学证明。不要给状态变量 `x`、`y` 赋数值后再把它们传给 `Compute`。

## 6. 查看结果和中间量

`FocusPeriod:-Show(ans)` 用于显示一份简洁结果。`ans` 本身是 Maple `table`，保留全部中间量，字符串键需要加双引号。

| Maple 访问方式 | 内容 |
|---|---|
| `ans["degree"]` | 输入多项式的最高总次数 |
| `ans["G"][j]`、`ans["H"][j]` | 极坐标系数，`ans["H"][0]=-1` |
| `ans["r_dot"]`、`ans["theta_dot"]` | 极坐标中的两个微分方程右端 |
| `ans["A"][j]` | 角速度倒数展开系数，`ans["A"][0]=-1` |
| `ans["R"][j]` | `dr/dtheta` 展开系数 |
| `ans["u"][n]`、`ans["du"][n]` | 径向解系数及其角度导数 |
| `ans["g"][n]`、`ans["p"][n]` | 焦点量、周期系数 |
| `ans["B"][n]` | 已处理方向符号的时间密度系数 |
| `ans["r_series"]`、`ans["P_series"]` | 截断径向解、截断周期展开 |
| `ans["time_density_series"]` | 代入径向解后的时间密度展开 |
| `ans["Ng"]`、`ans["Np"]` | 请求的焦点量和周期系数阶数 |
| `ans["u_order"]`、`ans["A_order"]` | 实际计算的中间阶数 |

只访问已经计算的阶数。例如 `Np=4` 时，不能把 `ans["p"][6]` 当作已计算结果。
`G` 的下标为 `0..degree`，`H` 的下标为 `0..degree-1`；线性系统的 `degree=1`，此时应查看 `ans["G"]`，而不要读取不存在的 `G[2]`。

角变量及半径变量使用模块私有符号，避免和系统中的同名参数碰撞。进行代入或提取系数时，请从结果表中获取这些符号：

```maple
th := ans["theta"]:
rh := ans["rho"]:

ans["G"][2];
ans["H"][0];
ans["A"][1];
ans["u"][2];
ans["g"][3];
ans["p"][4];

# 检查 u[2](0)=0。
simplify(eval(ans["u"][2], th=0));

# 从截断周期表达式中提取 rho^4 的系数。
coeff(ans["P_series"], rh, 4);
```

如果暂时只想核对极坐标化结果，可以跳过后续计算：

```maple
polar := FocusPeriod:-PolarData(X, Y, x, y):
polar["G"][2];
polar["H"][0];
polar["r_dot"];
polar["theta_dot"];
```

## 7. 高阶计算前需要注意什么？

符号表达式的大小可能随着阶数、参数数量和输入次数迅速增长。建议从 `Ng=3,Np=2` 或 `Ng=5,Np=4` 开始，确认输入和低阶结果后再提高阶数。

有参数关系时，先将已知关系代入系统；只需要具体参数结果时，先代入精确数值。打印大量高阶表达式也会占用时间，可先用冒号结束赋值语句，再单独查看需要的 `g[n]` 或 `p[n]`。

主程序采用截断多项式运算，只生成所需的幂次，但不保证任意高阶多参数系统都能在短时间内完成。若 Maple 无法求出某一步积分，程序会明确报错并停止，不把未求出的积分当作已计算系数。

程序中的 `g[n]` 保留完整返回映射系数，不自动施加前面各阶焦点量为零的条件。因此高阶结果可能含低阶系数的乘积；例如焦点示例的 `g[5]=a*b+3*Pi*a^2`。它们并非已经按低阶消失条件约化的 Lyapunov 常数。

程序也不自动删除奇数周期系数。中心示例的 `p[3]` 通常非零，是实际返回的结果。

## 8. 常见问题

| 现象或报错 | 原因与处理 |
|---|---|
| `read` 找不到文件 | 确认 `currentdir()` 显示的是项目根目录，目录内有 `src` 和 `examples`；检查拼写及 ZIP 是否解压 |
| 加载主程序后没有计算结果 | `src/focus_period.mpl` 只定义模块；继续调用 `Compute`，或执行一个示例文件 |
| 看不到赋值结果 | Maple 中冒号 `:` 抑制显示；用分号 `;` 显示表达式，或调用 `Show` |
| 提示线性矩阵必须为 `[[0,1],[-1,0]]` | 检查一次项是否恰为 `X=y+...`、`Y=-x+...`；程序不会自动翻转向量场或变换线性部分 |
| 提示原点必须是平衡点 | 输入中有非零常数项；应先自行确定平衡点并平移、整理为要求的标准系统 |
| 提示输入必须是多项式 | 不能直接输入含 `sin(x)`、`1/(1-x)`、时滞或微分方程等式的系统；当前程序处理的是多项式右端 |
| 提示状态变量类型错误或已有赋值 | 先执行 `restart`，重新加载主程序并定义系统；确保 `x`、`y` 为两个不同的未赋值变量 |
| 提示阶数不符合要求 | `Ng` 为 `0` 或不小于 `2` 的整数；`Np` 为非负整数，不能为小数、负数或 `Ng=1` |
| 无法给某个名字赋值 | 可能使用了 Maple 受保护名字；使用 `ans`、`th`、`rh` 等普通名字，保留 `Pi`、`I`、`D`、`O` 及内置函数名原来的用途 |
| 参数意外变成数值 | 工作表前面已经给该参数赋值；执行 `restart` 后从头运行，或用明确的 `eval` 参数集合代入 |
| 出现未能求积分的报错 | 先降低阶数、减少符号参数或代入已知关系，确认系统正确后再重试 |
| 升级文件后仍然显示旧版本 | Maple 内存中仍保存旧模块；执行 `restart` 后重新 `read "src/focus_period.mpl"` |

`restart` 会清除当前会话的变量赋值及已加载模块。执行后需要重新加载主程序、重新定义 `X` 和 `Y`。若修改了程序文件，仅保存磁盘文件不会更新已加载的模块。

路径建议用 Maple 字符串中的正斜杠，例如 `"D:/Maple/Maple_Focus_Period"`。避免把系统文件管理器中的反斜杠路径原样拼入多层字符串。

## 9. 文件结构与验证

```text
Maple_Focus_Period/
├─ README.md
├─ VERSION
├─ CHANGELOG.md
├─ AGENTS.md
├─ src/
│  └─ focus_period.mpl
├─ examples/
│  ├─ 01_center.mpl
│  ├─ 02_focus.mpl
│  ├─ user_system.mpl
│  └─ FocusPeriod_Worksheet.mw
├─ tests/
│  ├─ verify_maple.mpl
│  └─ verify_math.py
├─ tools/
│  └─ build_worksheet.mpl
└─ docs/
   ├─ mathematics.md
   └─ validation.md
```

普通使用者只需 Maple；Python 不是运行主程序的依赖。`verify_math.py` 是另用 SymPy 对数学递推与解析示例作独立核对的验证脚本。

在 Maple 中复核测试时，先将当前目录设为项目根目录，再执行：

```maple
restart:
currentdir("D:/Maple/Maple_Focus_Period"):
read "tests/verify_maple.mpl";
```

成功时最后会显示 `PASS` 及检查数量。验证覆盖解析中心、非中心径向示例、论文的完整 `p2` 公式、计算阶数选择、线性及高次输入、符号隔离和错误输入拒绝。实际运行环境和结果见 [docs/validation.md](docs/validation.md)。

## 10. 版本和后续更新

查看当前加载程序的版本：

```maple
FocusPeriod:-Version();
# v1.0.0 对应返回字符串 "1.0.0"。
```

仓库中的 [VERSION](VERSION) 保存版本号，[CHANGELOG.md](CHANGELOG.md) 记录每个版本的具体改动、接口变化和验证情况。固定版本标签便于今后复现同一组结果。

更新时，建议下载新的完整版本到单独文件夹，先保留自己的系统文件和旧版本结果。将旧版 `user_system.mpl` 中的 `X`、`Y`、`Ng`、`Np` 复制到新版入口，再按新版 README 设置路径、执行 `restart` 并重新加载。确认示例及自己的低阶结果后，再进行高阶计算。

本项目后续收到明确的修改任务时，会同步维护相关程序、说明、验证记录和版本变更。
每次完成维护后更新版本号；修复和说明更新递增修订版本，兼容功能扩展递增次版本，不兼容变更递增主版本。
具体维护约定保存在 [AGENTS.md](AGENTS.md)。本项目没有建立后台监测；提出修改需求后再执行相应更新。

个人模型和计算结果可放在自己创建的 `models/`、`results/` 文件夹，仓库默认忽略这些目录。
用 Git 管理项目时，拉取更新会保留这些未提交的个人文件；使用 ZIP 更新时，请自行保留并复制个人文件，不要直接覆盖整个旧目录。

## 11. 参考文献与数学约定

W. Zhang, X. Hou, and Z. Zeng, *Weak centers and bifurcation of critical periods in reversible cubic systems*, **Computers & Mathematics with Applications 40** (2000), 771–782。算法参考第 2 节式 (6)–(15) 及第 781–782 页附录。

原论文的标准系统线性部分为 `[[0,-1],[1,0]]`，对应 `H0=1`；本程序使用带 `H0` 的统一公式，保持使用者要求的 `H0=-1`。核对论文公式时，仅在验证示例中将论文整个向量场反号；这保留中心的闭轨及其正周期。主程序不会自动反号任何输入。

本仓库提供参考文献信息和公式核对说明，不包含论文 PDF。完整递推公式、方向解释及中心和形式周期的区别见 [docs/mathematics.md](docs/mathematics.md)。
