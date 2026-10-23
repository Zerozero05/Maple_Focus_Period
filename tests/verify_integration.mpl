# 公开接口组合、跨模块加载和参数分式的集成回归。
restart:
read "src/focus_period.mpl":
integration_checks := 0:
IntegrationCheck := proc(actual,expected,label::string)
    global integration_checks;
    if not evalb(actual=expected) then
        if not type(actual,algebraic) or not type(expected,algebraic)
           or not evalb(simplify(expand(actual-expected),trig)=0) then
            error "INTEGRATION FAILED: %1; actual=%2 expected=%3",label,actual,expected;
        end if;
    end if;
    integration_checks := integration_checks+1;
end proc:

Xc := (1+b*x)*(y+a*x*y): Yc := (1+b*x)*(-x+a*y^2):
before_reload := FocusPeriod:-Compute(Xc,Yc,x,y,5,4):
old_keys := ["degree","G","H","H0","r_dot","theta_dot","theta","r","rho",
             "Ng","Np","u_order","A_order","A","R","u","du","g","p","B",
             "r_series","P_series","time_density_series"]:
for key in old_keys do
    IntegrationCheck(assigned(before_reload[key]),true,cat("legacy field ",key)):
end do:

# 重新加载后仍能读取旧结果的私有符号；不得把新模块theta硬代入旧u。
read "src/focus_period.mpl":
after_reload := FocusPeriod:-Compute(Xc,Yc,x,y,5,4):
IntegrationCheck(eval(before_reload["u"][2],before_reload["theta"]=Pi),-2*a,"old angle survives reload"):
half_old := FocusPeriod:-HalfReturn(before_reload):
IntegrationCheck(half_old["composition_residual"],0,"new analysis accepts old module instance"):
period_old := FocusPeriod:-CenterPeriod(before_reload):
IntegrationCheck(period_old["tau"][1],Pi*b*(b-a),"amplitude analysis after module reload"):
IntegrationCheck(period_old["odd_residual"],0,"old result retains amplitude parity"):
scaling_old := FocusPeriod:-ScalingCheck(before_reload,2):
IntegrationCheck(scaling_old["status"],"verified","scaling after module reload"):
IntegrationCheck(FocusPeriod:-RunInfo(before_reload)["program_version"],"1.1.0","stored version after reload"):

# 参数分式是合法的多项式系数；共同因子检测须先处理分母。
rational_X := y+x^2/b+x*y:
rational_Y := -x+x^2+y^2/(b+1):
rational_direct := FocusPeriod:-Compute(rational_X,rational_Y,x,y,3,2):
rational_reference := FocusPeriod:-Compute(rational_X,rational_Y,x,y,3,2,
                                         method="quotient",factor_orbits=false):
IntegrationCheck(rational_direct["g"][3],rational_reference["g"][3],"rational parameter focus"):
IntegrationCheck(rational_direct["p"][2],rational_reference["p"][2],"rational parameter period"):
IntegrationCheck(rational_direct["input"]["exactness"],"exact","rational parameter exactness"):
rational_domain := FocusPeriod:-AnalyzeConditions(rational_direct,[b=0],[b]):
IntegrationCheck(rational_domain["p"][0]["domain_status"],"inconsistent","input domain retained for constant p0"):

function_coeff := FocusPeriod:-Compute(y+sin(a)*x*(x^2+y^2),
                                     -x+sin(a)*y*(x^2+y^2),x,y,3,0):
IntegrationCheck(function_coeff["g"][3],-sin(a),"nonpolynomial parameter coefficient remains constant in states"):
IntegrationCheck(before_reload["u_order"],5,"default radial range unchanged"):
IntegrationCheck(before_reload["A_order"],4,"default reciprocal range unchanged"):
IntegrationCheck(FocusPeriod:-RunInfo(before_reload)["factor_orbits"],true,"run records factor detection option"):
# 非中心形式时间转换实际含奇次项时，不能只对偶次投影求临界候选。
mixed_formal := FocusPeriod:-Compute(y+x^2+x*y+2*y^2+x^3,-x+x*y+y^2+x*y^2,x,y,5,4):
mixed_period := FocusPeriod:-CenterPeriod(mixed_formal,formal=true):
IntegrationCheck(mixed_period["odd_residual"],Pi*mixed_period["alpha"]^3,"noncenter odd amplitude term retained"):
odd_rejected := false:
try FocusPeriod:-CriticalPeriods(mixed_period): catch: odd_rejected := true: end try:
IntegrationCheck(odd_rejected,true,"full odd derivative cannot be discarded"):
printf("PASS: %d public interface integration checks.\n",integration_checks):
