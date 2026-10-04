# ======================================================================
# 自己的模型：先另存一份，例如 models/my_model.mpl，再修改 X、Y、Ng、Np。
# 从仓库根目录运行：read "examples/user_system.mpl";
# 先设置实际根目录：currentdir("D:/Maple/Maple_Focus_Period"):
# Windows 路径字符串使用 / 或 \\；单反斜杠会被当作转义符，导致找不到文件。
# 在可执行 Maple Input 中运行；一维代码无需转换成二维公式。
# 本模板的默认输入为线性中心，可立即运行；请替换为你的原始系统。
# ======================================================================
restart:
read "src/focus_period.mpl":

# 状态变量 x、y 与需要保留的符号参数应未赋值。
# 输入原始方程右端；原点为平衡点，线性部分必须为 [[0,1],[-1,0]]。
X := y:
Y := -x:

# Ng 是 rho 幂的最高指数；0 为关闭，或取 >=2；Np 可取 >=0。
Ng := 5:
Np := 4:

ans := FocusPeriod:-Compute(X,Y,x,y,Ng,Np):
FocusPeriod:-Show(ans);
FocusPeriod:-Version();

# 按需查询：不要对全局 theta、rho 直接代入私有符号。
th := ans["theta"]:
rh := ans["rho"]:
ans["G"];
ans["H"][0];
ans["A"][1];
ans["u"][2];
ans["g"][3];
ans["p"][2];
ans["P_series"];
eval(ans["u"][2],th=0);
