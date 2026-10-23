# 从仓库根目录执行：cmaple -q tests/verify_deep_core.mpl
# 原算法与新算法对照、独立隐式方程残差、次数支撑、轨道因子和局部参数环。
restart:
read "src/focus_period.mpl":
checks:=0:
Check:=proc(a,b,label::string)
global checks;
if not evalb(simplify(expand(a-b),trig)=0) then
error "FAILED %1 actual=%2 expected=%3",label,a,b;
end if;
checks:=checks+1;
end proc:
Reject:=proc(action::procedure,label::string)
global checks;
local rejected;
rejected:=false;
try action() catch: rejected:=true end try;
if not rejected then error "FAILED accepted %1",label end if;
checks:=checks+1;
end proc:
Compare:=proc(a::table,b::table,label::string)
local n;
for n from 2 to min(a["u_order"],b["u_order"]) do
Check(a["u"][n],b["u"][n],cat(label," u",n));
Check(a["du"][n],b["du"][n],cat(label," du",n));
end do;
for n from 2 to a["Ng"] do Check(a["g"][n],b["g"][n],cat(label," g",n)) end do;
for n from 0 to a["Np"] do
Check(a["p"][n],b["p"][n],cat(label," p",n));
Check(a["B"][n],b["B"][n],cat(label," B",n));
end do;
end proc:
Run:=proc()
global checks;
local X,Y,dc,qc,nf,sd,dd,cp,full,jet_result,ref,th,rh,n,jet0,jet2,
      rational_jet,linear,mixed,rr,res,namesafe;
X:=(1+b*x)*(y+a*x*y);Y:=(1+b*x)*(-x+a*y^2);
dc:=FocusPeriod:-Compute(X,Y,x,y,5,4);
qc:=FocusPeriod:-Compute(X,Y,x,y,5,4,method="quotient");
nf:=FocusPeriod:-Compute(X,Y,x,y,5,4,factor_orbits=false);
Compare(dc,qc,"center direct/quotient"); Compare(dc,nf,"center factor on/off");
Check(dc["orbit_factor"],1+b*x,"normalized factor");
Check(dc["orbit_degree"],2,"reduced orbit degree");
th:=dc["theta"];rh:=dc["rho"];
Check(dc["orbit_G"][2],a*sin(th),"orbit G2");
Check(dc["R"][3],0,"original quotient cancellation");
X:=y+a*x*(x^2+y^2)+b*y*(x^2+y^2);
Y:=-x+a*y*(x^2+y^2)-b*x*(x^2+y^2);
dc:=FocusPeriod:-Compute(X,Y,x,y,7,6);
qc:=FocusPeriod:-Compute(X,Y,x,y,7,6,method="quotient");
dd:=FocusPeriod:-Compute(X,Y,x,y,7,6,optimize=false);
Compare(dc,qc,"focus direct/quotient"); Compare(dc,dd,"focus sparse/dense");
# 分离变量积分给出 1/(2*a*r^2)-(b/a)*ln(r) 随 theta 线性增长；展开核对。
Check(dc["g"][7],-a*b^2-8*Pi*a^2*b-10*Pi^2*a^3,"focus exact g7");
cp:=FocusPeriod:-Compute(X,Y,x,y,0,6,compact=true);
full:=FocusPeriod:-Compute(X,Y,x,y,0,6);
Compare(cp,full,"compact period");
Check(cp["u_order"],5,"compact cubic u order");
Check(full["u_order"],6,"default retains old array extent");
X:=y+c*x^2*(x^2+y^2);Y:=-x+c*x*y*(x^2+y^2);
sd:=FocusPeriod:-Compute(X,Y,x,y,7,6);
dd:=FocusPeriod:-Compute(X,Y,x,y,7,6,optimize=false);
Compare(sd,dd,"quartic support");
Check(sd["u"][4],-c*sin(th),"quartic exact u4");
Check(sd["u"][7],2*c^2*sin(th)^2,"quartic exact u7");
if sd["orbit_structure"]["u_possible_orders"]<>[1,4,7] then
error "FAILED quartic structural support";
end if; checks:=checks+1;
X:=y+(a+eps)*x*(x^2+y^2)+(b+eps^2)*y*(x^2+y^2);
Y:=-x+(a+eps)*y*(x^2+y^2)-(b+eps^2)*x*(x^2+y^2);
jet_result:=FocusPeriod:-Compute(X,Y,x,y,7,6,jet=[eps,2]);
ref:=FocusPeriod:-Compute(X,Y,x,y,7,6);
for n from 2 to 7 do
Check(jet_result["u"][n],convert(series(ref["u"][n],eps=0,3),polynom),cat("jet u",n));
Check(jet_result["g"][n],convert(series(ref["g"][n],eps=0,3),polynom),cat("jet g",n));
if degree(jet_result["u"][n],eps)>2 then error "FAILED untruncated jet" end if;
end do;
for n from 0 to 6 do
Check(jet_result["p"][n],convert(series(ref["p"][n],eps=0,3),polynom),cat("jet p",n));
end do;
jet2:=FocusPeriod:-Compute(X,Y,x,y,5,4,method="quotient",jet=[eps,2]);
jet_result:=FocusPeriod:-Compute(X,Y,x,y,5,4,jet=[eps,2]);
Compare(jet_result,jet2,"jet direct/quotient");
jet0:=FocusPeriod:-Compute(X,Y,x,y,5,4,jet=[eps,0]);
ref:=FocusPeriod:-Compute(eval(X,eps=0),eval(Y,eps=0),x,y,5,4);
Compare(jet0,ref,"jet zero order");
rational_jet:=FocusPeriod:-Compute(y+x*(x^2+y^2)/(1+eps),-x+y*(x^2+y^2)/(1+eps),x,y,3,0,jet=[eps,2]);
Check(rational_jet["g"][3],-1+eps-eps^2,"analytic rational coefficient jet");
linear:=FocusPeriod:-Compute(y,-x,x,y,0,8,compact=true);
Check(linear["u_order"],1,"compact constant angular speed");
Check(linear["P_series"],2*Pi,"linear compact period");
# 混合二/三次、无公共因子、G/H 同时非零的非可分离模型。
X:=y+x^2+x*y+2*y^2+x^3;Y:=-x+x*y+y^2+x*y^2;
mixed:=FocusPeriod:-Compute(X,Y,x,y,7,6);
ref:=FocusPeriod:-Compute(X,Y,x,y,7,6,method="quotient",optimize=false,factor_orbits=false);
Compare(mixed,ref,"mixed generic direct/quotient");
Check(mixed["orbit_factor"],1,"mixed field has no common factor");
rr:=mixed["r"];
res:=convert(series(eval(mixed["theta_dot"],rr=mixed["r_series"])*
      diff(mixed["r_series"],th)-eval(mixed["r_dot"],rr=mixed["r_series"]),rh=0,8),polynom);
for n from 0 to 7 do Check(coeff(res,rh,n),0,cat("implicit radial residual",n)) end do;
res:=convert(series(eval(mixed["theta_dot"],rr=mixed["r_series"])*
      mixed["time_density_series"]+1,rh=0,7),polynom);
for n from 0 to 6 do Check(coeff(res,rh,n),0,cat("implicit time residual",n)) end do;
namesafe:=FocusPeriod:-Compute(y+(r+rho)*x*(x^2+y^2),-x+(r+rho)*y*(x^2+y^2),x,y,3,0);
Check(namesafe["g"][3],-r-rho,"parameters r and rho remain independent");
Reject(proc() FocusPeriod:-Compute(y,-x,x,y,3,2,method="unknown") end proc,"unknown method");
Reject(proc() FocusPeriod:-Compute(y,-x,x,y,3,2,jet=[x,2]) end proc,"state jet variable");
Reject(proc() FocusPeriod:-Compute(y,-x,x,y,3,2,jet=[eps,-1]) end proc,"negative jet order");
Reject(proc() FocusPeriod:-Compute(y+x^2/eps,-x,x,y,3,2,jet=[eps,2]) end proc,"singular jet coefficient");
Reject(proc() FocusPeriod:-PolarData(y+eps^2,-x,x,y,jet=[eps,1]) end proc,"jet cannot hide a shifted equilibrium");
Reject(proc() FocusPeriod:-PolarData((1+eps^2)*y,-x,x,y,jet=[eps,1]) end proc,"jet cannot hide a wrong linear part");
Reject(proc() FocusPeriod:-Compute(y+eps^2,-x,x,y,3,2,jet=[eps,1]) end proc,"Compute checks the original equilibrium");
jet_result:=FocusPeriod:-Compute(y+eps^2*x^3,-x,x,y,3,2,jet=[eps,1]);
Check(jet_result["degree"],1,"jet forwarded to polar extraction");
Check(jet_result["input"]["input_degree"],3,"original degree retained in metadata");
printf("PASS: %d deep core checks.\n",checks);
end proc:
Run();
