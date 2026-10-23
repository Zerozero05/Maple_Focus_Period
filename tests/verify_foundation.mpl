# 在仓库根目录运行：read "tests/verify_foundation.mpl":
read "src/focus_period.mpl":
# 先读取组装后的 FocusPeriod 模块，再读取本文件；本文件不改模块。
foundation_checks := 0:
FoundationAssert := proc(actual, expected, label::string)
    global foundation_checks;
    if not evalb(actual=expected) then
        error "FOUNDATION FAIL: %1; actual=%2 expected=%3",label,actual,expected;
    end if;
    foundation_checks := foundation_checks+1;
    return NULL;
end proc:
FoundationZero := proc(actual,label::string)
    FoundationAssert(simplify(expand(actual),trig),0,label);
end proc:
FoundationReject := proc(thunk::procedure,label::string)
    local rejected;
    rejected := false;
    try thunk(); catch: rejected := true; end try;
    FoundationAssert(rejected,true,label);
end proc:

foundation_input := FocusPeriod:-InputInfo(y+a*x^2/(b-1),-x,x,y):
FoundationAssert(foundation_input["exactness"],"exact","exact input"):
FoundationAssert(foundation_input["x"],x,"original x"):
FoundationAssert(foundation_input["y"],y,"original y"):
FoundationAssert({op(foundation_input["parameters"])},{a,b},"input coefficient names"):
FoundationAssert(foundation_input["denominator_nonzero"],[b-1],"visible denominator"):
FoundationAssert(FocusPeriod:-InputInfo(y+0.1*x^2,-x,x,y)["exactness"],"floating","floating coefficients"):
foundation_input := FocusPeriod:-InputInfo(v+theta*u^2+rho*u*v+r*v^2,-u,u,v):
FoundationAssert({op(foundation_input["parameters"])},{theta,rho,r},"input theta r rho remain parameters"):

FoundationAssert(FocusPeriod:-ZeroStatus(a-a)["status"],"zero","identity zero"):
FoundationAssert(FocusPeriod:-ZeroStatus(2*Pi)["status"],"nonzero","exact nonzero constant"):
FoundationAssert(FocusPeriod:-ZeroStatus(a)["status"],"unknown","generic nonzero polynomial is not pointwise nonzero"):
FoundationAssert(FocusPeriod:-ZeroStatus(a,nonzero=[a])["status"],"nonzero","explicit nonzero assumption"):
FoundationAssert(FocusPeriod:-ZeroStatus(a*(b-a),relations=[a*(b-a)=0],parameters=[a,b])["status"],"zero","ideal member"):
FoundationAssert(FocusPeriod:-ZeroStatus(b-a,relations=[a*(b-a)=0],parameters=[a,b])["status"],"unknown","do not lose a=0 branch"):
FoundationAssert(FocusPeriod:-ZeroStatus(b-a,relations=[a*(b-a)=0],parameters=[a,b],nonzero=[a])["status"],"unknown","saturation not silently assumed"):
FoundationAssert(FocusPeriod:-ZeroStatus(Pi*b*(b-a),relations=[b-a=0],parameters=[a,b])["status"],"zero","Pi is an exact coefficient"):
zs := FocusPeriod:-ZeroStatus(a/b,relations=[a=0],parameters=[a,b]):
FoundationAssert(zs["status"],"zero","zero rational numerator"):
FoundationAssert(zs["denominator_nonzero"],[b],"rational domain retained"):
zs := FocusPeriod:-ZeroStatus(a/b,relations=[b=0],parameters=[a,b]):
FoundationAssert(zs["status"],"unknown","undefined constrained quotient is not zero"):
FoundationAssert(zs["domain_status"],"inconsistent","denominator contradiction"):
zs := FocusPeriod:-ZeroStatus(a,relations=[a=0,a=1],parameters=[a]):
FoundationAssert(zs["domain_status"],"inconsistent","empty relation domain"):
FoundationAssert(zs["status"],"unknown","no vacuous center result on empty domain"):
FoundationAssert(FocusPeriod:-ZeroStatus(a,relations=[a^2=0],parameters=[a])["status"],"unknown","radical membership is not claimed"):
FoundationAssert(FocusPeriod:-ZeroStatus(0.1*a)["status"],"unknown","floating algebra not certified"):
FoundationReject(()->FocusPeriod:-ZeroStatus(a,relations=[a*b=0],parameters=[a]),"missing relation parameter rejected"):
FoundationReject(()->FocusPeriod:-ZeroStatus(a,parameters=[a,a]),"duplicate polynomial variables rejected"):

fc := FocusPeriod:-Compute((1+b*x)*(y+a*x*y),(1+b*x)*(-x+a*y^2),x,y,3,3):
fc_before := fc["p"][2]:
fa := FocusPeriod:-AnalyzeConditions(fc,[b-a=0],[a,b]):
FoundationAssert(fa["p"][2]["status"],"zero","condition analysis period reduction"):
FoundationZero(fc["p"][2]-fc_before,"raw p unchanged"):
FoundationAssert(fa["raw_result_unchanged"],true,"raw unchanged flag"):

fh := FocusPeriod:-Compute(y+a*x^2,-x-2*a*x*y,x,y,2,0):
cc := FocusPeriod:-CenterCertificate(fh):
FoundationAssert(cc["status"],"proved","Hamiltonian sufficient condition"):
FoundationAssert(cc["method"],"Hamiltonian","Hamiltonian method identified"):
FoundationAssert(cc["real_coefficients"],true,"real coefficient verification recorded"):
FoundationAssert(cc["assumptions"],[a::real],"real parameter premise explicit"):
cc := FocusPeriod:-CenterCertificate(fc):
FoundationAssert(cc["status"],"proved","reversible sufficient condition"):
FoundationAssert(cc["method"],"reversible_y_reflection","reversibility identified"):
ff := FocusPeriod:-Compute(y+a*x*(x^2+y^2),-x+a*y*(x^2+y^2),x,y,2,0):
FoundationAssert(ff["g"][2],0,"focus has zero finite low-order list"):
FoundationAssert(FocusPeriod:-CenterCertificate(ff)["status"],"unknown","finite zero focus quantities are not a proof"):
fz := FocusPeriod:-Compute(y+I*x^2,-x-2*I*x*y,x,y,2,0):
FoundationAssert(FocusPeriod:-CenterCertificate(fz)["status"],"unknown","complex Hamiltonian input is not a real center certificate"):
fz := FocusPeriod:-Compute(y+sqrt(a)*x^2,-x-2*sqrt(a)*x*y,x,y,2,0):
FoundationAssert(FocusPeriod:-CenterCertificate(fz)["status"],"unknown","unknown real coefficient domain not assumed away"):
fr := FocusPeriod:-Compute((1+x+y)*y,-(1+x+y)*x,x,y,2,0):
cc := FocusPeriod:-CenterCertificate(fr,first_integral=x^2+y^2):
FoundationAssert(cc["status"],"proved","verified polynomial first integral"):
FoundationAssert(cc["method"],"first_integral_extremum","first integral method identified"):
FoundationAssert(FocusPeriod:-CenterCertificate(fr,first_integral=-x^2-y^2)["status"],"proved","negative definite Hessian"):
FoundationAssert(FocusPeriod:-CenterCertificate(fr,first_integral=x^2-y^2)["status"],"unknown","invalid first integral rejected"):
FoundationAssert(FocusPeriod:-CenterCertificate(fr,first_integral=(x^2+y^2)^2)["status"],"unknown","degenerate Hessian not sufficient"):
FoundationAssert(FocusPeriod:-CenterCertificate(fr,first_integral=log(1+x^2+y^2))["status"],"unknown","unverified analytic expression not claimed"):
fz := FocusPeriod:-Compute(y+0.1*x^2,-x-0.2*x*y,x,y,2,0):
FoundationAssert(FocusPeriod:-CenterCertificate(fz)["status"],"unknown","floating center identity is not exact certificate"):

ri := FocusPeriod:-RunInfo(fc):
FoundationAssert(ri["source_identity_status"],"not_recorded","no invented source identity"):
FoundationAssert(ri["source_digest"],"","empty source digest"):
FoundationAssert(ri["input"]["X"],(1+b*x)*(y+a*x*y),"run retains original input"):
ri := FocusPeriod:-RunInfo(fc,source_file="src/focus_period.mpl",revision="local candidate"):
FoundationAssert(ri["source_digest_algorithm"],"sha1","actual supported digest algorithm"):
FoundationAssert(length(ri["source_digest"]),40,"SHA1 digest length"):
FoundationAssert(ri["source_digest"],FileTools:-Hash("src/focus_period.mpl",method="sha1"),"actual source bytes hashed"):
FoundationAssert(ri["source_identity_status"],"specified_file_hashed","hash scope is user-specified file"):
FoundationAssert(ri["revision"],"local candidate","user revision retained"):
FoundationReject(()->FocusPeriod:-RunInfo(fc,source_file="__missing_focus_period_source__.mpl"),"missing source rejected"):

fj := FocusPeriod:-ParameterJet(y+a*x*(x^2+y^2)+b*y*(x^2+y^2),-x+a*y*(x^2+y^2)-b*x*(x^2+y^2),x,y,5,4,{a=1,b=2},{a=1},1):
FoundationZero(fj["parameter_jet"]["g_coefficients"][3,0]+1,"jet focus base"):
FoundationZero(fj["parameter_jet"]["g_coefficients"][3,1]+1,"jet focus first derivative"):
FoundationZero(fj["parameter_jet"]["g_coefficients"][5,0]-(2+3*Pi),"jet nonlinear focus base"):
FoundationZero(fj["parameter_jet"]["g_coefficients"][5,1]-(2+6*Pi),"jet nonlinear focus derivative"):
FoundationZero(fj["parameter_jet"]["p_coefficients"][4,1]-8*Pi^2,"jet period derivative"):
FoundationAssert(degree(fj["g"][5],fj["epsilon"]),1,"jet is truncated during computation"):
fj0 := FocusPeriod:-ParameterJet(y+a*x*(x^2+y^2),-x+a*y*(x^2+y^2),x,y,3,0,{a=2},{a=1},0):
FoundationZero(fj0["g"][3]+2,"zero-order jet"):
FoundationAssert(fj0["parameter_jet"]["order"],0,"zero-order metadata"):
FoundationAssert(fj0["parameter_jet"]["original_input"]["X"],y+a*x*(x^2+y^2),"jet original input retained"):
fj2 := FocusPeriod:-ParameterJet(y+a*x*(x^2+y^2),-x+a*y*(x^2+y^2),x,y,5,0,{a=0},{a=1},2):
FoundationZero(fj2["parameter_jet"]["g_coefficients"][5,2]-3*Pi,"second-order jet coefficient"):
fjc := FocusPeriod:-ParameterJet(y+epsilon*x^2+a*x*(x^2+y^2),-x+a*y*(x^2+y^2),x,y,3,0,{a=0},{a=1},1):
FoundationAssert(evalb(fjc["epsilon"]=epsilon),false,"epsilon is private and does not collide with input parameter"):
FoundationReject(()->FocusPeriod:-ParameterJet(y,-x,x,y,2,0,{}, {a=1},1),"direction without base rejected"):
FoundationReject(()->FocusPeriod:-ParameterJet(y,-x,x,y,2,0,{x=0},{x=1},1),"state direction rejected"):
FoundationReject(()->FocusPeriod:-ParameterJet(y,-x,x,y,2,0,{a=x},{a=1},1),"state-dependent base rejected"):
FoundationReject(()->FocusPeriod:-ParameterJet(y,-x,x,y,2,0,{a=0,a=1},{a=1},1),"duplicate base rejected"):

# 原输入分母不能因输出恒为0而消失；运行身份使用计算当时的环境。
fdomain := FocusPeriod:-Compute((1+x/b)*y,-(1+x/b)*x,x,y,2,0):
fa := FocusPeriod:-AnalyzeConditions(fdomain,[b=0],[b]):
FoundationAssert(fa["g"][2]["domain_status"],"inconsistent","input denominator retained even for zero g"):
saved_digits := Digits:
Digits := saved_digits+3:
ri := FocusPeriod:-RunInfo(fc):
FoundationAssert(ri["Digits"],saved_digits,"record computation Digits rather than later inspection Digits"):
Digits := saved_digits:
FoundationAssert(ri["method"],"direct","actual engine recorded"):
fsaddle := table(["input"=FocusPeriod:-InputInfo(x,-y,x,y)]):
FoundationAssert(FocusPeriod:-CenterCertificate(fsaddle)["status"],"unknown","zero-divergence saddle is not an elliptic center"):
fshift := table(["input"=FocusPeriod:-InputInfo(1+y,-x,x,y)]):
FoundationAssert(FocusPeriod:-CenterCertificate(fshift)["status"],"unknown","shifted equilibrium cannot receive a center certificate"):
assume(assumed_parameter>0):
assumed_input := FocusPeriod:-InputInfo(y+assumed_parameter*x^2,-x,x,y):
FoundationAssert(nops(assumed_input["active_assumptions"]),1,"active worksheet assumptions are recorded"):
printf("FOUNDATION PASS: %d checks\n",foundation_checks):
