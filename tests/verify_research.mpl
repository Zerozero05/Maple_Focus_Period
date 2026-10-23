# 在仓库根目录执行；仅依赖单个核心 .mpl。
restart:
read "src/focus_period.mpl":
# 根目录执行，或先由集成脚本 read 带数学过程的模块。
checks_math := 0:
CheckMath := proc(actual,expected,label::string)
    global checks_math;
    if not evalb(simplify(expand(actual-expected),trig)=0) then
        error "FAILED math analysis: %1; actual=%2 expected=%3",label,actual,expected;
    end if;
    checks_math := checks_math+1;
end proc:
AssertMath := proc(condition::truefalse,label::string)
    global checks_math;
    if not condition then error "FAILED math analysis: %1",label; end if;
    checks_math := checks_math+1;
end proc:
RejectMath := proc(action::procedure,label::string)
    local rejected;
    rejected := false;
    try action(); catch: rejected := true; end try;
    AssertMath(rejected,label);
end proc:

RunMathAnalysisTests := proc()
    local Xc,Yc,center,half,forward,cp,critical,rh,al,zz,n,focus,linear,
          stable,unstable,formal_data,even_center,manual,impossible,low,
          input_data,realcritical,th,exact_period,mixed,mixed_half,compact,
          integral_center,jet_center,jet_focus,eps,jet_cp;
    Xc := (1+b*x)*(y+a*x*y);
    Yc := (1+b*x)*(-x+a*y^2);
    center := FocusPeriod:-Compute(Xc,Yc,x,y,5,4);

    rh := center["rho"];
    half := FocusPeriod:-HalfReturn(center);
    CheckMath(coeff(half["S_series"],rh,2),-2*a,"half-map s2");
    CheckMath(half["composition_residual"],0,"signed half-map composition");
    CheckMath(half["dependencies"]["g4_minus_s2_g3_over_2"],0,"g4 dependency");
    forward := FocusPeriod:-ForwardReturn(center);
    CheckMath(forward["return_series"],rh,"center forward return");
    CheckMath(forward["inverse_residual"],0,"center inverse residual");
    AssertMath(evalb(forward["stability"]="unknown"),"finite zero focus data do not prove stability");
    cp := FocusPeriod:-CenterPeriod(center);
    al := cp["alpha"];
    CheckMath(cp["alpha_series"],rh-a*rh^2+2*a^2*rh^3-4*a^3*rh^4,"natural amplitude");
    CheckMath(cp["inverse_residual"],0,"natural amplitude inverse");
    CheckMath(cp["odd_residual"],0,"center period has no odd natural amplitude terms");
    CheckMath(cp["tau"][1],Pi*b*(b-a),"first natural period constant");
    CheckMath(cp["tau"][2],Pi*b*(b-a)*(3*b^2-6*a*b-5*a^2)/4,"second natural period constant");
    AssertMath(evalb(cp["center_status"]="proved"),"independent reversible center proof");
    integral_center := FocusPeriod:-Compute((1+x+y)*y,-(1+x+y)*x,x,y,3,2);
    RejectMath(proc() FocusPeriod:-CenterPeriod(integral_center) end proc,"unrecognized center requires supplied evidence");
    jet_cp := FocusPeriod:-CenterPeriod(integral_center,first_integral=(x^2+y^2)/2);
    AssertMath(evalb(jet_cp["center_certificate"]["method"]="first_integral_extremum"),"first integral keyword reaches certificate");
    critical := FocusPeriod:-CriticalPeriods(cp,[a,b]);
    zz := critical["z"];
    CheckMath(critical["critical_polynomial"],cp["tau"][1]+2*cp["tau"][2]*zz,"critical polynomial");
    CheckMath(critical["jacobian"][1,1],-Pi*b,"period parameter Jacobian");
    CheckMath(critical["generic_rank"],2,"period generic rank");
    CheckMath(LinearAlgebra:-Rank(eval(critical["jacobian"],{a=0,b=0})),0,"generic rank may drop on a special branch");

    focus := FocusPeriod:-Compute(y+a*x*(x^2+y^2)+b*y*(x^2+y^2),
                                  -x+a*y*(x^2+y^2)-b*x*(x^2+y^2),x,y,5,4);
    forward := FocusPeriod:-ForwardReturn(focus);
    CheckMath(forward["g"][3],a,"forward g3 sign");
    CheckMath(forward["g"][5],3*Pi*a^2-a*b,"forward g5 nonlinear inversion");
    CheckMath(forward["inverse_residual"],0,"focus forward inverse");
    half := FocusPeriod:-HalfReturn(focus);
    CheckMath(half["composition_residual"],0,"focus half-map identity");
    AssertMath(evalb(forward["stability"]="unknown"),"parameter-dependent focus has unknown sign");
    mixed := FocusPeriod:-Compute(y+x*y+x*(x^2+y^2),-x+y^2+y*(x^2+y^2),x,y,4,0);
    mixed_half := FocusPeriod:-HalfReturn(mixed);
    CheckMath(mixed["g"][3],-1,"mixed quadratic cubic focus");
    CheckMath(mixed["g"][4],1,"nonzero even full return coefficient");
    CheckMath(mixed_half["dependencies"]["g4_minus_s2_g3_over_2"],0,"nontrivial g4 dependence");
    RejectMath(proc() FocusPeriod:-CenterPeriod(focus) end proc,"focus requires certificate");
    formal_data := FocusPeriod:-CenterPeriod(focus,formal=true);
    AssertMath(evalb(formal_data["center_status"]="formal_only"),"formal mode label");
    unstable := FocusPeriod:-ForwardReturn(FocusPeriod:-Compute(y+x*(x^2+y^2),-x+y*(x^2+y^2),x,y,3,0));
    stable := FocusPeriod:-ForwardReturn(FocusPeriod:-Compute(y-x*(x^2+y^2),-x-y*(x^2+y^2),x,y,3,0));
    AssertMath(evalb(unstable["stability"]="unstable"),"positive physical radial drift unstable");
    AssertMath(evalb(stable["stability"]="stable"),"negative physical radial drift stable");
    jet_center := FocusPeriod:-ParameterJet(Xc,Yc,x,y,5,4,{a=1,b=2},{a=1},1);
    eps := jet_center["epsilon"];
    jet_cp := FocusPeriod:-CenterPeriod(jet_center);
    CheckMath(jet_cp["inverse_residual"],0,"parameter jet amplitude inverse");
    CheckMath(jet_cp["odd_residual"],0,"parameter jet natural period parity");
    CheckMath(jet_cp["half_return"]["composition_residual"],0,"parameter jet half-return identity");
    CheckMath(jet_cp["tau"][1],2*Pi-2*Pi*eps,"parameter jet first period constant");
    CheckMath(jet_cp["tau"][2],-5*Pi/2-17*Pi*eps/2,"parameter jet second period constant");
    critical := FocusPeriod:-CriticalPeriods(jet_cp);
    AssertMath(evalb(critical["jet"]=jet_center["jet"]),"critical candidate preserves parameter truncation metadata");
    jet_focus := FocusPeriod:-ParameterJet(y+a*x*(x^2+y^2),-x+a*y*(x^2+y^2),x,y,5,0,{a=0},{a=1},1);
    eps := jet_focus["epsilon"];
    forward := FocusPeriod:-ForwardReturn(jet_focus);
    CheckMath(forward["g"][3],eps,"parameter jet forward g3");
    CheckMath(forward["g"][5],0,"inverse must not recreate omitted parameter powers");
    CheckMath(forward["inverse_residual"],0,"parameter jet forward inverse residual");

    # 独立解析周期核对到六阶；不从同一递推表达式反推期望值。
    even_center := FocusPeriod:-Compute((1+2*x)*(y+x*y),(1+2*x)*(-x+y^2),x,y,6,6);
    cp := FocusPeriod:-CenterPeriod(even_center);
    al := cp["alpha"];
    exact_period := convert(series(4*Pi*sqrt((1+sqrt(1+4*al^2))/2)-2*Pi,al=0,7),polynom);
    CheckMath(cp["P_alpha_series"],exact_period,"independent exact period in natural amplitude");
    CheckMath(cp["odd_residual"],0,"sixth order center natural parity");
    CheckMath(cp["tau"][3],21*Pi/4,"sixth order natural period constant");

    linear := FocusPeriod:-Compute(y,-x,x,y,0,0);
    cp := FocusPeriod:-CenterPeriod(linear);
    CheckMath(cp["P_alpha_series"],2*Pi,"zero order center period");
    critical := FocusPeriod:-CriticalPeriods(cp);
    CheckMath(critical["critical_polynomial"],0,"zero order no critical information");
    CheckMath(critical["generic_rank"],0,"empty parameter rank");
    AssertMath(evalb(critical["degree_status"]="no_information_beyond_computed_order"),"constant truncation is not isochrony proof");
    RejectMath(proc() FocusPeriod:-CenterPeriod(center,5) end proc,"uncomputed period order rejected");
    low := copy(center);
    low["u_order"] := 2;
    RejectMath(proc() FocusPeriod:-CenterPeriod(low) end proc,"insufficient radial order rejected");
    compact := FocusPeriod:-Compute(y*(1+x^2+y^2),-x*(1+x^2+y^2),x,y,0,4,compact=true);
    AssertMath(evalb(compact["u_order"]=3),"compact period uses fewer radial orders");
    RejectMath(proc() FocusPeriod:-CenterPeriod(compact) end proc,"actual compact result cannot silently supply missing radial terms");
    RejectMath(proc() FocusPeriod:-CriticalPeriods(cp,[a,a]) end proc,"repeated Jacobian parameter rejected");
    RejectMath(proc() FocusPeriod:-CriticalPeriods(table()) end proc,"missing center-period metadata rejected");

    # 独立构造已知截断多项式，明确区分正、负候选和截断含义。
    manual := table(["tau_order"=2,"tau"=table([1=-2,2=1]),"odd_residual"=0,"center_status"="proved"]);
    critical := FocusPeriod:-CriticalPeriods(manual);
    AssertMath(evalb(critical["candidates"][1]["positivity"]="positive"),"positive critical candidate");
    CheckMath(critical["candidates"][1]["z"],1,"critical candidate z");
    manual["tau"][1] := 2;
    critical := FocusPeriod:-CriticalPeriods(manual);
    AssertMath(evalb(critical["candidates"][1]["positivity"]="nonpositive"),"negative candidate excluded");
    manual["tau"][1] := -2*a;
    critical := FocusPeriod:-CriticalPeriods(manual);
    AssertMath(evalb(critical["candidates"][1]["positivity"]="undetermined"),"parameter-dependent positivity stays unknown");
    manual["odd_residual"] := a;
    RejectMath(proc() FocusPeriod:-CriticalPeriods(manual) end proc,"odd amplitude terms cannot be silently discarded");
    printf("PASS: %d mathematical analysis checks.\n",checks_math);
end proc:
RunMathAnalysisTests():
