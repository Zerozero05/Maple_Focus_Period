# 从仓库根目录运行：read "examples/02_focus.mpl";
restart:
read "src/focus_period.mpl":

# 径向三次焦点：检验方向符号及未经低阶消失条件约化的焦点量。
# 此时 P_series 仅表示用户约定下的形式时间展开。
X := y+a*x*(x^2+y^2)+b*y*(x^2+y^2):
Y := -x+a*y*(x^2+y^2)-b*x*(x^2+y^2):
Ng := 5:
Np := 4:
ans := FocusPeriod:-Compute(X,Y,x,y,Ng,Np):
FocusPeriod:-Show(ans);

# 预期：G[2]=0，G[3]=a，H[1]=0，H[2]=-b。
# g[2]=g[4]=0，g[3]=-a，g[5]=a*b+3*Pi*a^2。
# p[1]=p[3]=0，p[2]=-2*Pi*b，p[4]=2*Pi*b^2+4*Pi^2*a*b。
