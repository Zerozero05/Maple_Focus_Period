# 从仓库根目录运行：read "examples/01_center.mpl";
# 先设置实际根目录：currentdir("D:/Maple/Maple_Focus_Period"):
# Windows 路径字符串使用 / 或 \\；单反斜杠会被当作转义符，导致找不到文件。
# 在可执行 Maple Input 中运行；一维代码无需转换成二维公式。
# 独立示例会 restart 清除工作表变量。核心算法在 src/focus_period.mpl。
restart:
read "src/focus_period.mpl":

# 非等时解析中心：精确解 r=rho/(1+a*rho*(1-cos(theta)))。
# 这个中心示例同时有非零 G[2]、G[3]、H[1]，适合核对周期混合项。
X := (1+b*x)*(y+a*x*y):
Y := (1+b*x)*(-x+a*y^2):
Ng := 5:
Np := 4:
ans := FocusPeriod:-Compute(X,Y,x,y,Ng,Np):
FocusPeriod:-Show(ans);

# 预期：G[2]=a*sin(theta)，G[3]=a*b*cos(theta)*sin(theta)，
# H[0]=-1，H[1]=-b*cos(theta)，H[2]=0，g[2]..g[5]=0。
# p[1]=0；p[2]=Pi*b*(b-a)；p[3]=-2*Pi*a*b*(b-a)；
# p[4]=3*Pi*b*(b-a)*(b^2-2*a*b+5*a^2)/4。
# Maple 可能显示三角恒等或展开后的等价表达式。

th := ans["theta"]:
rh := ans["rho"]:
ans["A"][0];
ans["A"][1];
ans["u"][2];
eval(ans["u"][2],th=0);
eval(ans["P_series"],{a=1,b=2});
# 数值参数代入后：2*Pi+2*Pi*rh^2-4*Pi*rh^3+15*Pi*rh^4/2。
