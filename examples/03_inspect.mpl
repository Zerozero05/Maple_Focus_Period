# ======================================================================
# 查看关键变量：在仓库根目录运行 read "examples/03_inspect.mpl";
# 先用 currentdir("D:/Maple/Maple_Focus_Period"): 设置实际项目根目录。
# Windows 路径使用 / 或双反斜杠；红色一维 Maple Input 可以直接执行。
# 本示例会 restart；完整含义与下标范围见 README 第6节。
# ======================================================================
restart:
read "src/focus_period.mpl":

# 已知中心示例，含二次和三次项，适合查看全部递推量。
X := (1+b*x)*(y+a*x*y):
Y := (1+b*x)*(-x+a*y^2):
Ng := 5:
Np := 4:
ans := FocusPeriod:-Compute(X,Y,x,y,Ng,Np):
FocusPeriod:-Show(ans);

# 输入次数 d、实际径向最高阶 Nu、倒数最高阶 NA。
d := ans["degree"]:
Nu := ans["u_order"]:
NA := ans["A_order"]:
[d, ans["Ng"], ans["Np"], Nu, NA];  # [3,5,4,5,4]

# 查看完整对象；结果会保留原有 Array 下标和 table 结构。
ans["G"];
ans["H"];
ans["A"];
ans["R"];
ans["u"];
ans["du"];
eval(ans["g"]); # g 为 table，用 eval 显示实际表项。
ans["p"];
ans["B"];

# 明确按下标顺序显示为列表。G[0]=G[1]=0，R[0]=R[1]=0。
[seq(ans["G"][j], j=0..d)];
[seq(ans["H"][j], j=0..d-1)];
[seq(ans["A"][j], j=0..NA)];
[seq(ans["R"][j], j=0..Nu)];
[seq(ans["u"][n], n=1..Nu)];
[seq(ans["du"][n], n=1..Nu)];
[seq(ans["g"][n], n=2..ans["Ng"])]; # Ng=0 时为空列表。
[seq(ans["p"][n], n=0..ans["Np"])];
[seq(ans["B"][n], n=0..ans["Np"])];

# 查看极坐标方程和已经截断的展开式（表达式本身不附 O 项）。
ans["r_dot"];
ans["theta_dot"];
ans["r_series"];
ans["time_density_series"];
ans["P_series"];

# 从结果中获取真正使用的私有符号，再代入角度、半径或参数。
th := ans["theta"]:
rr := ans["r"]:
rh := ans["rho"]:
ans["G"][2];                          # a*sin(th)
simplify(eval(ans["G"][2], th=Pi/2)); # a
simplify(eval(ans["r_dot"], {th=Pi/2, rr=1})); # a
simplify(eval(ans["u"][2], th=0));    # 0，不能写 u[2](0)。
simplify(eval(ans["u"][2], th=2*Pi)); # 本中心示例为 0。

# 核对同一系数的定义：返回位移 / (2*Pi)、时间密度积分。
simplify(eval(ans["u"][3], th=2*Pi)/(2*Pi)-ans["g"][3]); # 0
simplify(int(ans["B"][2], th=0..2*Pi)-ans["p"][2]);     # 0
simplify(ans["B"][2]+ans["A"][1]*ans["u"][2]+ans["A"][2]); # 0

# 参数集中代入，不改变 ans 里的原符号结果。
eval(ans["p"][2], {a=1,b=2});       # 2*Pi
eval(ans["P_series"], {a=1,b=2});   # 2*Pi+2*Pi*rh^2-4*Pi*rh^3+15*Pi*rh^4/2
coeff(ans["P_series"], rh, 4);      # 与 ans["p"][4] 相同。

# 只计算极坐标数据时使用 PolarData；它没有 A、u、g、p 等后续结果。
polar_data := FocusPeriod:-PolarData(X,Y,x,y):
polar_data["G"];
polar_data["H"];
polar_data["H0"];                  # -1
