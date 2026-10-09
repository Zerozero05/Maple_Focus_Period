# 验证说明

本记录对应程序版本 `v1.0.1` 和 **2026-10-10** 的说明更新。算法与已发布 v1.0.1 一致；变量查询和工作表已在 Maple 2025.1 验证。

## Maple 程序实际运行

环境：Maple 2025.1，Windows 64 位。

从仓库根目录运行：

```maple
restart:
currentdir("D:/Maple/Maple_Focus_Period"):
read "tests/verify_maple.mpl";
```

或在已安装 Maple 的终端中执行 `cmaple -q tests/verify_maple.mpl`。
若 `cmaple` 未配置到 PATH，使用你实际 Maple 安装目录下的内核完整路径。

预期结束信息：

```text
PASS: 58 Maple checks, including the published p2 formula.
```

核对范围包括：

- 固定角方向、`G[2]`、`G[3]`、`H[0]`、`H[1]`、`H[2]`。
- 解析中心精确径向解到 `rho^5`，积分初值和周期系数到 `p[4]`。
- 径向三次焦点的 `g[3]`、完整 `g[5]` 和含漂移项的形式 `p[4]`。
- Zhang 等（2000）第 774 页完整七参数 `p[2]` 公式。
- 两类输出独立指定、关闭其中一类、线性系统和四次系统。
- 状态变量换名、参数名与模块私有符号隔离。
- 错误线性部分、非零常数项、非多项式/有理输入、不合法阶数的拒绝。

`examples/01_center.mpl`、`examples/02_focus.mpl` 和 `examples/user_system.mpl`
均为可执行示例。v1.0.1 完善关键变量说明和查询示例、更新版本查询，没有改动递推算法。

## 关键变量查询

`examples/03_inspect.mpl` 在 Maple 2025.1 中完整执行成功。README 第 6 节的三组
查询代码也逐条执行通过：整体对象、按下标列表、角度与参数代入、初值与返回值、
截断表达式系数提取，以及只查看极坐标结果的入口。

另作 140 项查询和范围核对，覆盖返回键与类型、`Nu`/`NA`、`G`/`H`/`A`/`R`/
`u`/`du`/`g`/`B`/`p` 的下标和定义、`Ng=0`、`Np=0`、纯线性系统及符号隔离。

本次实际查询发现旧说明使用了受保护名字 `polar`；已改为 `polar_data`。
`g` 是 `table`，直接输入 `ans["g"]` 可能只显示内部表名，因此说明采用
`eval(ans["g"])` 或按已计算下标生成列表。`u[n]` 是表达式，初值通过 `eval` 查询。

特别注意：命令行 Maple 在报错后可能仍然返回退出码 0。因此验证时须同时检查
输出中没有 `Error`，并确实出现上述 PASS 信息。

## 独立数学核对

`tests/verify_math.py` 使用精确解和独立符号计算核对数学公式，不需要调用 Maple。
已核对 Python 3.12 / SymPy 1.14.0 环境。

```text
python -m pip install -r requirements-test.txt
python tests/verify_math.py
```

核对范围为中心和焦点的极坐标方程、倒数与径向递推、焦点量、周期混合项及论文
的完整七参数 `p[2]`。最后应出现：

```text
All independent symbolic checks passed.
```

GitHub Actions 每次推送或提出 Pull Request 时运行这一核对。GitHub 运行环境没有
安装商业 Maple 内核，所以绿色状态代表独立数学检查通过，不能替代 `.mpl` 文件
在真实 Maple 中的运行验证。

## 工作表模板

`examples/FocusPeriod_Worksheet.mw` 使用原生 Maple 工作表格式，调用
`src/focus_period.mpl`，没有在工作表内复制算法。其可执行输入可通过原生 Maple
文档工具读取、导出为 Maple 代码，再由内核执行核对。下载后须修改首个输入块中的
项目路径，然后按顺序运行。Windows 路径字符串推荐使用 `/`；使用反斜杠时，
每个分隔符必须写成 `\\`。红色一维输入是可执行 Maple Input。

本版本已在 Maple 2025.1 完成原生读取及整页执行，13 项断言通过，
核对版本显示、中心 `G[2]`、`G[3]`、`H[0]`、`H[1]`、`g[3]`、`p[2..4]`、
积分初值，以及焦点 `g[3]`、`g[5]` 和形式 `p[4]`。
维护工作表内容时，可在仓库根目录执行 `read "tools/build_worksheet.mpl":`
重新生成模板，再作同样的读取与运行检查。

## 路径和单文件入口

在真实 Maple 2025.1 中复现了单反斜杠被转义而导致目录不存在的问题；
已核对 README 中正斜杠和双反斜杠两种路径写法的解析结果一致。

仅复制 `focus_period.mpl` 到独立临时目录，在该目录中执行 `read` 并计算中心示例，
版本、`H[0]`、`g[3]`、`p[2]`、`p[4]` 共 5 项检查通过。单文件入口不依赖项目
其他文件；完整 `.mw` 模板则仍需按其相对路径保留 `src` 文件夹。

## 数学结论的边界

测试给出具体模型、公式和输入检查的证据，不构成任意阶算法的形式化证明。
有限阶 `g[n]=0` 也不足以证明中心。一般焦点下 `P_series` 是规定方向的形式时间
展开；在独立确认中心的条件下才把它解释为闭轨的正周期。
