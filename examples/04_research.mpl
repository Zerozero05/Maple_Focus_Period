# v1.1.0开发候选的研究分析示例；在仓库根目录执行。
# Windows路径必须用 / 或双反斜杠；一维Maple Input可直接执行。
restart:
read "src/focus_period.mpl":
X := (1+b*x)*(y+a*x*y):
Y := (1+b*x)*(-x+a*y^2):
ans := FocusPeriod:-Compute(X,Y,x,y,5,4):

# 原G/H/A/R不变；轨道因子只简化求u，不从时间密度中删除。
ans["orbit_factor"];
ans["structure"]["u_possible_orders"];
ans["input"]["exactness"];

# 当前可逆中心有独立充分条件，不依赖有限g消失。
certificate := FocusPeriod:-CenterCertificate(ans):
certificate["status"]; certificate["method"];
half_data := FocusPeriod:-HalfReturn(ans):
half_data["S_series"]; half_data["composition_residual"];
period_data := FocusPeriod:-CenterPeriod(ans):
period_data["alpha_series"];
period_data["tau"][1]; period_data["tau"][2];
period_data["odd_residual"];
critical_data := FocusPeriod:-CriticalPeriods(period_data,[a,b]):
critical_data["critical_polynomial"];
critical_data["jacobian"]; critical_data["generic_rank"];
# 候选根是截断结果；generic秩在特殊参数分支上可能下降。

conditions := FocusPeriod:-AnalyzeConditions(ans,[b-a=0],[a,b]):
conditions["p"][2]["status"];
# 原ans["p"][2]仍为Pi*b*(b-a)，不会被条件约化覆盖。
jet_data := FocusPeriod:-ParameterJet(X,Y,x,y,5,4,{a=1,b=2},{a=1,b=0},1):
jet_data["p"][2];
jet_data["parameter_jet"]["p_coefficients"][2,1];
# 一阶系数为-2*Pi；K阶Taylor系数是方向K阶导数/K!。

cartesian := FocusPeriod:-Lyapunov(X,Y,x,y,6):
eval(cartesian["L"]); cartesian["residual"];
scaling := FocusPeriod:-ScalingCheck(ans,2):
scaling["status"];
run_info := FocusPeriod:-RunInfo(ans,source_file="src/focus_period.mpl"):
run_info["source_digest_algorithm"]; run_info["source_digest"];

focus_data := FocusPeriod:-Compute(y+a*x*(x^2+y^2)+b*y*(x^2+y^2),
                 -x+a*y*(x^2+y^2)-b*x*(x^2+y^2),x,y,5,4):
forward_data := FocusPeriod:-ForwardReturn(focus_data):
forward_data["g"][3]; forward_data["g"][5];
# 分别为a和3*Pi*a^2-a*b；高阶系数不是全部逐项反号。
