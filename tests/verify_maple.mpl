# 在项目目录执行：cmaple -q tests/verify_maple.mpl
# 对解析解、论文公布的独立公式及输入边界作核对。
restart:
read "src/focus_period.mpl":

checks := 0:
Check := proc(actual, expected, label::string)
    global checks;
    if not evalb(simplify(expand(actual-expected),trig)=0) then
        error "FAILED: %1; actual=%2; expected=%3", label, actual, expected;
    end if;
    checks := checks+1;
end proc:

Reject := proc(action::procedure, label::string)
    global checks;
    local rejected;
    rejected := false;
    try
        action();
    catch:
        rejected := true;
    end try;
    if not rejected then
        error "FAILED: invalid input was accepted: %1", label;
    end if;
    checks := checks+1;
end proc:

RunTests := proc()
local Xc, Yc, center, th, rh, exact_radius, n, Xf, Yf, focus,
      Xp, Yp, paper, paper_p2, only_g, only_p, linear_result, none,
      fourth, renamed, namesafe;
global checks;

# 1. 解析中心：与精确径向解展开比较，非零二次项及奇数周期项。
Xc := (1+b*x)*(y+a*x*y):
Yc := (1+b*x)*(-x+a*y^2):
center := FocusPeriod:-Compute(Xc,Yc,x,y,5,4):
th := center["theta"]:
rh := center["rho"]:
Check(center["H0"],-1,"fixed angular sign"):
Check(center["G"][2],a*sin(th),"center G2"):
Check(center["G"][3],a*b*cos(th)*sin(th),"center G3"):
Check(center["H"][1],-b*cos(th),"center H1"):
Check(center["H"][2],0,"center H2"):
exact_radius := convert(series(rh/(1+a*rh*(1-cos(th))),rh=0,6),polynom):
Check(center["r_series"],exact_radius,"exact center radial solution"):
for n from 2 to 5 do
    Check(center["g"][n],0,cat("center g",n)):
    Check(eval(center["u"][n],th=0),0,cat("initial condition u",n)):
end do:
Check(center["p"][0],2*Pi,"positive period constant"):
Check(center["p"][1],0,"center p1"):
Check(center["p"][2],Pi*b*(b-a),"center p2"):
Check(center["p"][3],-2*Pi*a*b*(b-a),"center p3"):
Check(center["p"][4],3*Pi*b*(b-a)*(b^2-2*a*b+5*a^2)/4,"center p4"):

# 2. 径向三次焦点：g 的符号及含漂移项的形式时间系数。
Xf := y+a*x*(x^2+y^2)+b*y*(x^2+y^2):
Yf := -x+a*y*(x^2+y^2)-b*x*(x^2+y^2):
focus := FocusPeriod:-Compute(Xf,Yf,x,y,5,4):
Check(focus["G"][3],a,"focus G3"):
Check(focus["H"][2],-b,"focus H2"):
Check(focus["u"][3],-a*th,"focus u3"):
Check(focus["u"][5],a*b*th+3*a^2*th^2/2,"focus u5"):
Check(focus["g"][3],-a,"focus g3 orientation"):
Check(focus["g"][5],a*b+3*Pi*a^2,"focus g5 without reduction"):
Check(focus["p"][2],-2*Pi*b,"formal time p2"):
Check(focus["p"][4],2*Pi*b^2+4*Pi^2*a*b,"formal time p4 drift term"):

# 3. Zhang等(2000)第774页的完整七参数 p2，整体反号以满足用户标准。
# 时间反向不改变闭轨的正周期；没有更改极坐标。
Xp := y-a1*x^2-a2*y^2-a3*x^2*y-a4*y^3:
Yp := -x-b1*x*y-b2*x^3-b3*x*y^2:
paper := FocusPeriod:-Compute(Xp,Yp,x,y,3,2):
Check(paper["G"][2],-a1*cos(th)^3-(a2+b1)*sin(th)^2*cos(th),"paper G2"):
Check(paper["G"][3],-(a3+b2)*cos(th)^3*sin(th)-(a4+b3)*cos(th)*sin(th)^3,"paper G3"):
Check(paper["H"][1],(a1-b1)*cos(th)^2*sin(th)+a2*sin(th)^3,"paper H1"):
Check(paper["H"][2],(a3-b3)*cos(th)^2*sin(th)^2-b2*cos(th)^4+a4*sin(th)^4,"paper H2"):
paper_p2 := Pi*(-5*a1*b1+4*a1^2-3*b3+3*a3+b1^2-a2*b1+10*a1*a2-9*b2+9*a4+10*a2^2)/12:
Check(paper["p"][2],paper_p2,"published seven-parameter p2"):
Check(paper["g"][2],0,"reversible paper g2"):
Check(paper["g"][3],0,"reversible paper g3"):

# 4. 不同输出阶数、关闭某类输出、线性系统和四次系统。
only_g := FocusPeriod:-Compute(Xf,Yf,x,y,5,0):
only_p := FocusPeriod:-Compute(Xc,Yc,x,y,0,4):
Check(only_g["g"][5],focus["g"][5],"focus-only computation"):
Check(only_g["A_order"],3,"minimal focus reciprocal order"):
Check(only_g["P_series"],2*Pi,"period disabled"):
Check(only_p["p"][4],center["p"][4],"period-only computation"):
Check(nops([indices(only_p["g"])]),0,"no unrequested focus output"):
linear_result := FocusPeriod:-Compute(y,-x,x,y,4,3):
Check(linear_result["degree"],1,"linear degree"):
Check(linear_result["P_series"],2*Pi,"linear period"):
Check(linear_result["r_series"],rh,"linear radius"):
none := FocusPeriod:-Compute(y,-x,x,y,0,0):
Check(none["u_order"],1,"both outputs disabled"):
Check(none["A_order"],0,"zero reciprocal order"):
fourth := FocusPeriod:-Compute(y*(1+c*x^3),-x*(1+c*x^3),x,y,2,6):
Check(fourth["degree"],4,"degree is not fixed at three"):
Check(fourth["G"][4],0,"fourth degree circular center"):
Check(fourth["H"][3],-c*cos(th)^3,"fourth degree H3"):
Check(fourth["p"][6],5*Pi*c^2/8,"fourth degree exact circular period"):

# 状态变量换名，以及输入参数名恰好是 theta、r、rho 的情况。
renamed := FocusPeriod:-Compute(v,-z,z,v,2,2):
Check(renamed["P_series"],2*Pi,"renamed state variables"):
namesafe := FocusPeriod:-Compute(y+theta*x*(x^2+y^2),-x+theta*y*(x^2+y^2),x,y,3,0):
Check(namesafe["g"][3],-theta,"parameter theta is independent of private angle"):

# 5. 无效输入必须明确拒绝。
Reject(proc() FocusPeriod:-Compute(-y,x,x,y,3,2) end proc,"wrong rotation"):
Reject(proc() FocusPeriod:-Compute(1+y,-x,x,y,3,2) end proc,"shifted equilibrium"):
Reject(proc() FocusPeriod:-Compute(y+sin(x),-x,x,y,3,2) end proc,"nonpolynomial"):
Reject(proc() FocusPeriod:-Compute(y/(1-x),-x,x,y,3,2) end proc,"rational vector field"):
Reject(proc() FocusPeriod:-Compute(y,-x,x,y,1,2) end proc,"g1 ambiguity"):
Reject(proc() FocusPeriod:-Compute(y,-x,x,y,-1,2) end proc,"negative order"):
Reject(proc() FocusPeriod:-Compute(y,-x,x,y,3,3/2) end proc,"fractional order"):
Reject(proc() FocusPeriod:-Compute(y,-x,x,x,3,2) end proc,"same state variable"):

printf("PASS: %d Maple checks, including the published p2 formula.\n",checks);
end proc:
RunTests();
