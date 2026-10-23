# 独立笛卡尔同调内核及坐标缩放恒等式核对；在仓库根目录执行。
restart:
read "src/focus_period.mpl":

checks_cartesian := 0:
CheckCartesian := proc(actual,expected,label::string)
    global checks_cartesian;
    if not evalb(simplify(expand(actual-expected),trig)=0) then
        error "FAILED Cartesian verification: %1; actual=%2 expected=%3",label,actual,expected;
    end if;
    checks_cartesian := checks_cartesian+1;
end proc:
AssertCartesian := proc(condition::truefalse,label::string)
    global checks_cartesian;
    if not condition then error "FAILED Cartesian verification: %1",label; end if;
    checks_cartesian := checks_cartesian+1;
end proc:
RejectCartesian := proc(action::procedure,label::string)
    local rejected;
    rejected := false;
    try action(); catch: rejected := true; end try;
    AssertCartesian(rejected,label);
end proc:

RunCartesianTests := proc()
    local linear,radial,quintic,center,mixed,renamed,namesafe,
          radial_polar,quintic_polar,mixed_polar,Xm,Ym,m,n,
          raw_residual,degree_residual,component,scaled,center_polar,
          symbolic_scale,weight_test,jet_result,epsilon,jet_scaled;
    linear := FocusPeriod:-Lyapunov(y,-x,x,y,2);
    CheckCartesian(linear["V"],x^2+y^2,"linear minimal potential");
    CheckCartesian(linear["residual"],0,"linear minimal residual");
    AssertCartesian(evalb(nops([indices(linear["L"])])=0),"N2 has no requested Lyapunov quantities");
    linear := FocusPeriod:-Lyapunov(y,-x,x,y,9);
    CheckCartesian(linear["V"],x^2+y^2,"linear higher potential");
    for n from 1 to 3 do
        CheckCartesian(linear["L"][n],0,cat("linear L",n));
    end do;
    CheckCartesian(linear["residual"],0,"linear higher residual");

    # 解析径向系统：高阶完整g包含前面焦点量的乘积，不能充当高阶L。
    radial := FocusPeriod:-Lyapunov(y+a*x*(x^2+y^2)+b*y*(x^2+y^2),
                                    -x+a*y*(x^2+y^2)-b*x*(x^2+y^2),x,y,8);
    CheckCartesian(radial["L"][1],2*a,"radial cubic L1");
    CheckCartesian(radial["L"][2],0,"radial cubic L2");
    CheckCartesian(radial["L"][3],0,"radial cubic L3");
    CheckCartesian(radial["V"],x^2+y^2,"radial exact potential");
    radial_polar := FocusPeriod:-Compute(y+a*x*(x^2+y^2)+b*y*(x^2+y^2),
                                        -x+a*y*(x^2+y^2)-b*x*(x^2+y^2),x,y,5,0);
    CheckCartesian(radial_polar["g"][3],-radial["L"][1]/2,"first nonzero g3 and L1");
    CheckCartesian(radial_polar["g"][5],a*b+3*Pi*a^2,"complete g5 retains nonlinear lower-order contributions");
    AssertCartesian(not evalb(radial_polar["g"][5]=-radial["L"][2]/2),"raw g5 differs from independent L2");

    quintic := FocusPeriod:-Lyapunov(y+c*x*(x^2+y^2)^2,-x+c*y*(x^2+y^2)^2,x,y,8);
    CheckCartesian(quintic["L"][1],0,"quintic lower L1 vanishes");
    CheckCartesian(quintic["L"][2],2*c,"quintic first L2");
    CheckCartesian(quintic["L"][3],0,"quintic L3");
    quintic_polar := FocusPeriod:-Compute(y+c*x*(x^2+y^2)^2,-x+c*y*(x^2+y^2)^2,x,y,5,0);
    CheckCartesian(quintic_polar["g"][5],-c,"quintic g5 orientation");
    CheckCartesian(quintic_polar["g"][5],-quintic["L"][2]/2,"first nonzero g5 and L2");

    center := FocusPeriod:-Lyapunov((1+b*x)*(y+a*x*y),(1+b*x)*(-x+a*y^2),x,y,8);
    for n from 1 to 3 do
        CheckCartesian(center["L"][n],0,cat("reversible center L",n));
    end do;
    CheckCartesian(center["residual"],0,"reversible center homological residual");

    # 一般混合系统：直接对返回V代回原场取齐次分量，独立复核内置residual。
    Xm := y+x^2+2*x*y+y^2+x^3+x*y^2;
    Ym := -x+x*y+2*y^2-x^2*y;
    mixed := FocusPeriod:-Lyapunov(Xm,Ym,x,y,8);
    raw_residual := expand(diff(mixed["V"],x)*Xm+diff(mixed["V"],y)*Ym-mixed["target_series"]);
    degree_residual := 0;
    for m from 0 to 8 do
        degree_residual := degree_residual
            +add(coeff(coeff(raw_residual,x,m-n),y,n)*x^(m-n)*y^n,n=0..m);
    end do;
    CheckCartesian(degree_residual,0,"independent total-degree homological residual");
    CheckCartesian(mixed["residual"],0,"reported mixed homological residual");
    for m from 4 by 2 to 8 do
        component := mixed["V_components"][m];
        CheckCartesian(coeff(component,x,m),0,cat("even-degree normalization m",m));
    end do;
    mixed_polar := FocusPeriod:-Compute(Xm,Ym,x,y,3,2);
    CheckCartesian(mixed_polar["g"][3],-mixed["L"][1]/2,"mixed first g3-L1 relation");

    renamed := FocusPeriod:-Lyapunov(v+a*z*(z^2+v^2),-z+a*v*(z^2+v^2),z,v,6);
    CheckCartesian(renamed["L"][1],2*a,"renamed states preserve L1");
    CheckCartesian(renamed["residual"],0,"renamed states residual");
    namesafe := FocusPeriod:-Lyapunov(y+(theta+r+rho)*x*(x^2+y^2),
                                      -x+(theta+r+rho)*y*(x^2+y^2),x,y,6);
    CheckCartesian(namesafe["L"][1],2*(theta+r+rho),"parameters theta r rho remain independent");
    CheckCartesian(namesafe["L"][2],0,"namesafe independent higher L");

    RejectCartesian(proc() FocusPeriod:-Lyapunov(y,-x,x,y,1) end proc,"N below two");
    RejectCartesian(proc() FocusPeriod:-Lyapunov(y,-x,x,y,0) end proc,"zero N");
    RejectCartesian(proc() FocusPeriod:-Lyapunov(-y,x,x,y,4) end proc,"opposite linear rotation");
    RejectCartesian(proc() FocusPeriod:-Lyapunov(y+1,-x,x,y,4) end proc,"shifted origin");
    RejectCartesian(proc() FocusPeriod:-Lyapunov(y+sin(x),-x,x,y,4) end proc,"nonpolynomial field");
    RejectCartesian(proc() FocusPeriod:-Lyapunov(y,-x,x,x,4) end proc,"same state names");

    center_polar := FocusPeriod:-Compute((1+b*x)*(y+a*x*y),(1+b*x)*(-x+a*y^2),x,y,5,4);
    scaled := FocusPeriod:-ScalingCheck(center_polar,2);
    AssertCartesian(evalb(scaled["status"]="verified"),"center scaling verified");
    for n from 2 to 5 do
        CheckCartesian(scaled["g_residual"][n],0,cat("center scaling g",n));
    end do;
    for n from 0 to 4 do
        CheckCartesian(scaled["p_residual"][n],0,cat("center scaling p",n));
    end do;
    symbolic_scale := FocusPeriod:-ScalingCheck(center_polar,s);
    AssertCartesian(evalb(symbolic_scale["status"]="verified"),"symbolic coordinate scaling verified");
    AssertCartesian(evalb(symbolic_scale["scale_nonzero"]=(s<>0)),"symbolic scale nonzero condition retained");
    scaled := FocusPeriod:-ScalingCheck(mixed_polar,3);
    AssertCartesian(evalb(scaled["status"]="verified"),"mixed field scaling verified");
    CheckCartesian(scaled["g_residual"][3],0,"mixed g has weight n-1");
    CheckCartesian(scaled["p_residual"][2],0,"mixed p has weight n");
    weight_test := FocusPeriod:-Compute(y+a*x*(x^2+y^2),-x+a*y*(x^2+y^2),x,y,5,0);
    scaled := FocusPeriod:-ScalingCheck(weight_test,s);
    CheckCartesian(scaled["g_residual"][5],0,"symbolic high focus weight");
    jet_result := FocusPeriod:-ParameterJet(y+a*x*(x^2+y^2),-x+a*y*(x^2+y^2),
                                           x,y,5,0,{a=1},{a=1},1);
    jet_scaled := FocusPeriod:-ScalingCheck(jet_result,2);
    AssertCartesian(evalb(jet_scaled["status"]="verified"),"parameter jet scaling verified in coefficient ring");
    CheckCartesian(jet_scaled["g_residual"][5],0,"jet scaling high-order focus weight");
    scaled := FocusPeriod:-ScalingCheck(weight_test,2.0);
    AssertCartesian(evalb(scaled["status"]="approximate_check"),"floating scale does not claim exact verification");
    RejectCartesian(proc() FocusPeriod:-ScalingCheck(center_polar,0) end proc,"zero scale");
    RejectCartesian(proc() FocusPeriod:-ScalingCheck(center_polar,x) end proc,"scale depends on state");
    RejectCartesian(proc() FocusPeriod:-ScalingCheck(center_polar,center_polar["theta"]) end proc,"scale depends on private angle");
    printf("PASS: %d Cartesian and scaling checks.\n",checks_cartesian);
end proc:
RunCartesianTests():
