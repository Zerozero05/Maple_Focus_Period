# ======================================================================
# 标准平面多项式系统：焦点量与周期系数
# v1.1.0；完整六参数调用与原始g/p定义保持兼容。
# 输入：x' = X(x,y), y' = Y(x,y)，线性部分必须是 [[0,1],[-1,0]]。
# 固定极坐标：x=r*cos(theta), y=r*sin(theta)，因此 H[0]=-1。
#
# 用法（本文件只定义模块，不执行 restart，也不覆盖用户的系统）：
#   currentdir("D:/Maple/Maple_Focus_Period"):
#   read "src/focus_period.mpl":
#   ans := FocusPeriod:-Compute(X, Y, x, y, 5, 4):
#   FocusPeriod:-Show(ans);
# Windows 路径：Maple 字符串中推荐使用 /，如 "D:/Maple/focus_period.mpl"。
# 若使用反斜杠，每个目录分隔符必须写成 \\，如 "D:\\Maple\\focus_period.mpl"。
# 不要直接粘贴资源管理器的单反斜杠路径；它会被当作转义符处理。
# currentdir 失败会使 read 失败，继而出现 FocusPeriod 不是模块的错误；先修正路径。
# 单文件下载后可直接 read "D:/Maple/focus_period.mpl":，无需其他项目文件。
# 红色的一维 Maple Input 是可执行代码，不必转换成二维公式；普通文字区不执行。
# 最后两个参数分别为焦点量最高指数 Ng、周期系数最高指数 Np。
# Ng=0 关闭焦点量输出；Ng>=2 输出 g[2]..g[Ng]；Np=0 只保留 P(0)。
# 六个位置参数保持兼容；可选关键字：method="direct"（默认）或 "quotient"，
# optimize=true（跳过可证明的零阶），factor_orbits=true（约去共同轨道因子），
# compact=false（默认保留旧中间量范围），jet=[] 或 [epsilon,K]（局部参数展开）。
# 改用 compact=true 后请按 ans["u_order"] 查询实际数组范围。
#
# 关键量都保存在 Compute 返回的 table 中，例如 ans["G"][2]：
# G[j] 是 dr/dt 的 r^j 系数，H[j] 是 dtheta/dt 的 r^j 系数。
# A[j] 是 1/(dtheta/dt) 的 r^j 系数；R[j] 是 dr/dtheta 的 r^j 系数。
# u[n] 是 r(theta,rho) 的 rho^n 系数；du[n] 是它对 theta 的导数。
# g[n] 是返回位移的归一化系数；B[n] 是带方向负号的时间密度的 rho^n 系数。
# p[n] 是 B[n] 在 0..2*Pi 的积分；不能直接用 -int(A[n],...) 替代。
# 查看一整组用 ans["A"];，逐项列表用 [seq(ans["A"][j],j=0..ans["A_order"])];
# u[n] 保存的是表达式，不是函数；初值用 eval(ans["u"][n],ans["theta"]=0)。
# 下标范围、所有结果字段和查询示例见 README 第6节、examples/03_inspect.mpl。
#
# 本程序采用用户约定：g[n]=u[n](2*Pi)/(2*Pi)，n>=2。
# u[1]=1 是初始半径项，不是返回位移，故不定义 g[1]。
# 所有 g[n] 均为完整返回系数，不自动模掉前面的焦点量。
#
# 周期约定：P(rho)=-int(1/H(theta,r(theta,rho)),theta=0..2*Pi)。
# 中心情形：这是同一闭轨的正周期，P(0)=2*Pi。
# 一般焦点情形：这里只是逆时间绕行一圈的形式时间展开，不是真实周期。
# 正时间对应 theta:0 -> -2*Pi；本程序不改变角变量，也不翻转 g[n]。
# 有限个 g[n]=0 不能单独证明中心，中心条件须另行分析。
#
# 参考：Zhang, Hou, Zeng (2000), Computers & Mathematics with
# Applications 40, 771-782，第2节式(6)--(15)，附录第781--782页。
# 原文 H[0]=1；此处推广到任意多项式次数并使用用户的 H[0]=-1。
# ======================================================================

FocusPeriod := module()
    option package;
    export Compute, PolarData, Show, Version, InputInfo, ZeroStatus,
           AnalyzeConditions, CenterCertificate, RunInfo, ParameterJet,
           Lyapunov, ScalingCheck, HalfReturn, ForwardReturn,
           CenterPeriod, CriticalPeriods;
    local Clean, Cut, JetCut, CheckJet, IntegralFromZero, Reciprocal,
          RadialCoefficients, RadialSolution, PeriodSeries, StructureData,
          OrbitData, theta, r, rho, VisibleDenominators, RelationPolynomial,
          ReduceCondition, HomogeneousPart, MathCut, MathCompose,
          MathInverse, MathConstantSign, alpha, z;

    # theta、r、rho 是模块私有符号，避免与输入中的同名参数混淆。
    # 需要代入时，请使用 ans["theta"]、ans["r"]、ans["rho"]。

    # 每个发布版本同步维护此版本号、根目录 VERSION 和 CHANGELOG.md。
    Version := proc()
        return "1.1.0";
    end proc;

    # jet=[] 保持普通符号计算；[epsilon,K] 表示在 epsilon=0 的 K 阶系数环中运算。
    # 截断发生在每次乘积、化简与积分前后，不是先生成完整参数通式再截断。
    CheckJet := proc(jet::list, x::name, y::name)
        if jet=[] then return NULL end if;
        if nops(jet)<>2 then
            error "jet must be [] or [epsilon, nonnegative integer]";
        end if;
        if not type(jet[1],name) or not type(jet[2],nonnegint) then
            error "jet must be [] or [epsilon, nonnegative integer]";
        end if;
        if member(jet[1],{x,y,theta,r,rho}) then
            error "The jet variable must be independent of the state and polar variables";
        end if;
        return NULL;
    end proc;

    JetCut := proc(f, jet::list:=[])
        local e, K, s, v, K1;
        if jet=[] or evalb(f=0) then return f end if;
        e := jet[1]; K := jet[2];
        # 分母在基点为零、负幂或不可展开分支均明确报错，不能丢弃奇异项。
        try
            s := series(f,e=0,K+1);
            v := convert(s,polynom);
        catch:
            error "Cannot form a regular parameter jet at %1=0: %2",e,f;
        end try;
        if not type(v,polynom(anything,e)) then
            error "The parameter expansion is not a regular Taylor jet: %1",f;
        end if;
        return add(coeff(v,e,K1)*e^K1,K1=0..K);
    end proc;

    Clean := proc(f, jet::list:=[])
        return JetCut(simplify(expand(JetCut(f,jet)),trig),jet);
    end proc;

    # 对 rho 多项式精确截断，避免直接展开高次幂产生不需要的项。
    Cut := proc(f, N::nonnegint, jet::list:=[])
        local e, k;
        e := expand(JetCut(f,jet));
        return add(JetCut(coeff(e,rho,k),jet)*rho^k,k=0..N);
    end proc;

    # 从 0 积分到 upper；减去积分常数，保证 u[n](0)=0。
    # 若 Maple 留下未计算的积分，明确停止，不将其冒充为已求出的系数。
    IntegralFromZero := proc(f, upper, jet::list:=[])
        local primitive, integrand;
        if evalb(f=0) then
            return 0;
        end if;
        integrand := JetCut(f,jet);
        primitive := int(integrand, theta);
        if has(primitive, {'int', 'Int'}) then
            error "Maple could not evaluate an integral; integrand: %1", f;
        end if;
        return Clean(eval(primitive, theta=upper)
                     - eval(primitive, theta=0),jet);
    end proc;

    # 第1步：检查原始系统，自动极坐标化并按 r 的幂提取 G、H。
    PolarData := proc(X, Y, x::name, y::name, {jet::list:=[]})
        local XX, YY, d, xp, yp, radial, angular, G, H, j;
        CheckJet(jet,x,y);
        if evalb(x=y) then
            error "The two state variables must be distinct";
        end if;
        if not type(X, polynom(anything, {x,y}))
           or not type(Y, polynom(anything, {x,y})) then
            error "X and Y must be polynomials in the state variables";
        end if;
        # 先核验完整输入；不能用参数截断掩盖高阶epsilon中的平衡点/线性错误。
        XX := expand(X);
        YY := expand(Y);
        if not evalb(Clean(eval(XX, {x=0,y=0}))=0)
           or not evalb(Clean(eval(YY, {x=0,y=0}))=0) then
            error "The origin must be an equilibrium";
        end if;
        if not evalb(Clean(eval(diff(XX,x), {x=0,y=0}))=0)
           or not evalb(Clean(eval(diff(XX,y), {x=0,y=0}))=1)
           or not evalb(Clean(eval(diff(YY,x), {x=0,y=0}))=-1)
           or not evalb(Clean(eval(diff(YY,y), {x=0,y=0}))=0) then
            error "The linear matrix must be [[0,1],[-1,0]]";
        end if;

        XX := expand(JetCut(XX,jet));
        YY := expand(JetCut(YY,jet));
        d := max(degree(XX, {x,y}), degree(YY, {x,y}));
        xp := subs({x=r*cos(theta), y=r*sin(theta)}, XX);
        yp := subs({x=r*cos(theta), y=r*sin(theta)}, YY);
        radial := expand(xp*cos(theta) + yp*sin(theta));
        angular := expand((yp*cos(theta) - xp*sin(theta))/r);

        G := Array(0..d, fill=0);
        H := Array(0..d-1, fill=0);
        for j from 0 to d do
            G[j] := Clean(coeff(radial, r, j),jet);
        end do;
        for j from 0 to d-1 do
            H[j] := Clean(coeff(angular, r, j),jet);
        end do;
        return table([
            "degree"=d, "G"=G, "H"=H, "H0"=H[0],
            "r_dot"=add(G[j]*r^j, j=2..d),
            "theta_dot"=add(H[j]*r^j, j=0..d-1),
            "theta"=theta, "r"=r, "rho"=rho
        ]);
    end proc;

    # 第2步：1/(H[0]+H[1]*r+...) = sum(A[n]*r^n,n>=0)。
    # A[0]=1/H[0]；A[n]=-sum(H[j]*A[n-j])/H[0]。
    Reciprocal := proc(H::Array, d::posint, N::nonnegint, jet::list:=[])
        local A, n, j;
        A := Array(0..N, fill=0);
        A[0] := 1/H[0];
        for n from 1 to N do
            A[n] := Clean(-add(H[j]*A[n-j], j=1..min(n,d-1))/H[0],jet);
        end do;
        return A;
    end proc;

    # 第3步：dr/dtheta=sum(R[j]*r^j,j>=2)，
    # R[j]=sum(G[k]*A[j-k],k=2..min(j,d))。
    # r(theta,rho)=rho+sum(u[n](theta)*rho^n,n>=2)。
    # 求 u[n]' 时只需要 u[1]..u[n-1]，因为径向方程从 r^2 开始。
    RadialCoefficients := proc(G::Array, A::Array, d::posint,
                               N::posint, jet::list:=[])
        local R, j, k;
        R := Array(0..N, fill=0);
        for j from 2 to N do
            R[j] := Clean(add(G[k]*A[j-k],k=2..min(j,d)),jet);
        end do;
        return R;
    end proc;

    # 次数支撑是必要条件：只跳过可由次数证明为零的阶，不用数值猜测。
    # 从 G 的至少一个非零阶出发，再加入 G/H 的次数移位，得到 u 的支撑上界。
    StructureData := proc(G::Array,H::Array,d::posint,Nu::posint,Np::nonnegint)
        local gs, hs, shifts, reachable, up, bp, M, n, j, q;
        gs := [seq(`if`(evalb(G[j]=0),NULL,j-1),j=2..d)];
        hs := [seq(`if`(evalb(H[j]=0),NULL,j),j=1..d-1)];
        shifts := {op(gs),op(hs)};
        M := max(Nu-1,Np);
        reachable := Array(0..M,fill=false); reachable[0] := true;
        for n from 1 to M do
            for j in shifts do
                if j<=n and reachable[n-j] then
                    reachable[n] := true; break;
                end if;
            end do;
        end do;
        up := Array(1..Nu,fill=false); up[1] := true;
        for n from 2 to Nu do
            for j in gs do
                if j<=n-1 and reachable[n-1-j] then
                    up[n] := true; break;
                end if;
            end do;
        end do;
        bp := Array(0..Np,fill=false); bp[0] := true;
        for n from 1 to Np do
            for j in hs do
                if j<=n and reachable[n-j] then
                    bp[n] := true; break;
                end if;
            end do;
        end do;
        q := 0; if nops(hs)>0 then q := min(op(hs)) end if;
        return table(["radial_shifts"=gs,"angular_shifts"=hs,
            "u_possible"=up,"period_possible"=bp,"first_angular_order"=q,
            "u_possible_orders"=[seq(`if`(up[n],n,NULL),n=1..Nu)],
            "period_possible_orders"=[seq(`if`(bp[n],n,NULL),n=0..Np)]]);
    end proc;

    # 同时约去 X/Y 的公共多项式因子。只简化轨道计算；原场及时间密度不变。
    # 因子归一化为 q(0)=1；除法后必须仍为状态变量多项式。
    OrbitData := proc(X,Y,x::name,y::name,enabled::truefalse,
                       original::table,jet::list:=[])
        local q, q0, XX, YY, data;
        XX := JetCut(X,jet); YY := JetCut(Y,jet); q := 1;
        if enabled then
            # 系数可含参数分式；先取规范化分子，不能把分式直接交给多项式gcd。
            q := gcd(numer(normal(XX)),numer(normal(YY)));
            q0 := normal(eval(q,{x=0,y=0}));
            if evalb(q0=0) then
                error "A common orbit factor vanishes at the origin";
            end if;
            q := normal(q/q0);
            XX := normal(XX/q); YY := normal(YY/q);
            if not type(XX,polynom(anything,{x,y}))
               or not type(YY,polynom(anything,{x,y})) then
                error "The common orbit factor did not yield polynomial vector fields";
            end if;
        end if;
        if evalb(q=1) then
            # 未约去因子时复用已核验的极坐标数据，避免重复展开原场。
            data := copy(original);
        else
            data := PolarData(XX,YY,x,y,':-jet'=jet);
        end if;
        data["factor"] := q;
        return data;
    end proc;

    RadialSolution := proc(G::Array,H::Array,A::Array,d::posint,N::posint,
                            method::string,optimize::truefalse,support::table,
                            jet::list:=[])
        local R, u, du, solution, power, rhs, n, j, k, powers, angular;
        u := Array(1..N, fill=0);
        du := Array(1..N, fill=0);
        u[1] := 1;
        solution := rho;
        R := Array(0..N,fill=0);
        if method="quotient" then R := RadialCoefficients(G,A,d,N,jet) end if;
        for n from 2 to N do
            if optimize and not support["u_possible"][n] then next end if;
            rhs := 0;
            if method="quotient" then
                power := solution;
                for j from 2 to n do
                    power := Cut(power*solution,n,jet);
                    rhs := JetCut(rhs+R[j]*coeff(power,rho,n),jet);
                end do;
            else
                # H(s)*s_theta=G(s)：当前 n 阶仅依赖已经求出的低阶解。
                powers := Array(0..min(d,n),fill=0); powers[0] := 1;
                for j from 1 to min(d,n) do
                    powers[j] := Cut(powers[j-1]*solution,n,jet);
                end do;
                for j from 2 to min(d,n) do
                    rhs := JetCut(rhs+G[j]*coeff(powers[j],rho,n),jet);
                end do;
                angular := 0;
                for j from 1 to min(d-1,n) do
                    angular := JetCut(angular+H[j]*coeff(
                        Cut(powers[j]*diff(solution,theta),n,jet),rho,n),jet);
                end do;
                rhs := JetCut((rhs-angular)/H[0],jet);
            end if;
            du[n] := Clean(rhs,jet);
            u[n] := IntegralFromZero(du[n],theta,jet);
            solution := solution + u[n]*rho^n;
        end do;
        return table(["R"=R, "u"=u, "du"=du, "r_series"=solution]);
    end proc;

    # 第4步：先将 r(theta,rho) 代入 sum(A[j]*r^j)，再提取 rho 系数。
    # 注意不能直接对 A[n] 积分；A[1]*u[n] 等混合项同样贡献 p[n]。
    # 负号只作用在时间密度上，保持原极坐标与 u[n] 完全不变。
    PeriodSeries := proc(A::Array,H::Array,d::posint,solution,N::nonnegint,
                         method::string,optimize::truefalse,support::table,
                         jet::list:=[])
        local p, B, rp, power, density, j, n, h;
        p := Array(0..N, fill=0);
        B := Array(0..N, fill=0);
        h := Array(0..N,fill=0);
        rp := Cut(solution,N,jet);
        power := 1;
        if method="quotient" then
            density := A[0];
            for j from 1 to N do
                power := Cut(power*rp,N,jet);
                density := JetCut(density+A[j]*power,jet);
            end do;
            density := -Cut(density,N,jet);
        else
            # 先复合原场角速度，再由 H(r(theta,rho))*B=-1 逐阶求倒数。
            density := H[0];
            for j from 1 to min(d-1,N) do
                power := Cut(power*rp,N,jet);
                density := JetCut(density+H[j]*power,jet);
            end do;
            density := Cut(density,N,jet);
            for n from 0 to N do h[n] := Clean(coeff(density,rho,n),jet) end do;
            B[0] := -1/H[0];
            for n from 1 to N do
                if optimize and not support["period_possible"][n] then next end if;
                B[n] := Clean(-add(h[j]*B[n-j],j=1..n)/H[0],jet);
            end do;
        end if;
        for n from 0 to N do
            if method="quotient" then B[n] := Clean(coeff(density,rho,n),jet) end if;
            p[n] := IntegralFromZero(B[n],2*Pi,jet);
        end do;
        return table([
            "p"=p, "B"=B, "composed_H"=h,
            "time_density_series"=add(B[n]*rho^n, n=0..N),
            "P_series"=add(p[n]*rho^n, n=0..N)
        ]);
    end proc;

    # 主接口：两个最高阶独立指定；这里只计算实际需要的中间阶数。
    # u 需要到 max(Ng,Np,1)，A 需要到 max(Np,N_u-2,0)。
    # 默认 direct 用交叉相乘径向递推；quotient 保留旧算法用于独立核对。
    # optimize=false 禁用次数跳过，factor_orbits=false 禁用公共因子约去。
    # 默认仍保留旧版中间量范围；compact=true 才缩减仅周期所需的 u 阶数。
    Compute := proc(X, Y, x::name, y::name,
                    Ng::nonnegint, Np::nonnegint,
                    {method::string:="direct",optimize::truefalse:=true,
                     factor_orbits::truefalse:=true,compact::truefalse:=false,
                     jet::list:=[]})
        local result, d, Nu, NA, A, radial, period, g, n, original_support,
              orbit, orbit_support, orbit_A, period_order, q, old_Nu;
        if Ng=1 then
            error "Use Ng=0 to omit focus quantities, or Ng>=2";
        end if;
        if not member(method,{"direct","quotient"}) then
            error "method must be direct or quotient";
        end if;
        CheckJet(jet,x,y);
        result := PolarData(X,Y,x,y,':-jet'=jet);
        d := result["degree"];
        old_Nu := max(1,Ng,Np);
        original_support := StructureData(result["G"],result["H"],d,old_Nu,Np);
        q := original_support["first_angular_order"];
        period_order := 1;
        if q>0 and Np>=q then period_order := max(1,Np-q+1) end if;
        Nu := old_Nu;
        if compact then Nu := max(1,Ng,period_order) end if;
        NA := max(0, Np, Nu-2);
        A := Reciprocal(result["H"],d,NA,jet);
        orbit := OrbitData(X,Y,x,y,factor_orbits,result,jet);
        orbit_support := StructureData(orbit["G"],orbit["H"],orbit["degree"],Nu,Np);
        # direct不需要轨道场的倒数展开；quotient且因子为1时复用原A。
        orbit_A := A;
        if method="quotient" and not evalb(orbit["factor"]=1) then
            orbit_A := Reciprocal(orbit["H"],orbit["degree"],max(0,Nu-2),jet);
        end if;
        radial := RadialSolution(orbit["G"],orbit["H"],orbit_A,
                    orbit["degree"],Nu,method,optimize,orbit_support,jet);
        period := PeriodSeries(A,result["H"],d,radial["r_series"],Np,
                               method,optimize,original_support,jet);
        g := table();
        for n from 2 to Ng do
            g[n] := Clean(eval(radial["u"][n],theta=2*Pi)/(2*Pi),jet);
        end do;

        result["Ng"] := Ng;
        result["Np"] := Np;
        result["u_order"] := Nu;
        result["A_order"] := NA;
        result["A"] := A;
        # 可查询的 A/R 永远对应用户原始场，不因轨道因子优化而改变定义。
        result["R"] := RadialCoefficients(result["G"],A,d,Nu,jet);
        result["u"] := radial["u"];
        result["du"] := radial["du"];
        result["g"] := g;
        result["p"] := period["p"];
        result["B"] := period["B"];
        result["r_series"] := radial["r_series"];
        result["P_series"] := period["P_series"];
        result["time_density_series"] := period["time_density_series"];
        result["method"] := method;
        result["optimize"] := optimize;
        result["factor_orbits"] := factor_orbits;
        result["compact"] := compact;
        result["jet"] := jet;
        result["requested_u_order"] := old_Nu;
        result["period_u_order"] := period_order;
        result["orbit_factor"] := orbit["factor"];
        result["orbit_G"] := orbit["G"];
        result["orbit_H"] := orbit["H"];
        result["orbit_degree"] := orbit["degree"];
        result["structure"] := original_support;
        result["orbit_structure"] := orbit_support;
        result["input"] := InputInfo(X,Y,x,y);
        if method="direct" then result["composed_H"] := period["composed_H"] end if;
        return result;
    end proc;

    # 简洁输出；更大的中间表达式保存在结果表中，按需自行查看。
    Show := proc(result::table)
        local j;
        printf("degree = %d, H0 = %a\n", result["degree"], result["H0"]);
        for j from 2 to result["degree"] do
            printf("G[%d] = %a\n", j, result["G"][j]);
        end do;
        for j from 1 to result["degree"]-1 do
            printf("H[%d] = %a\n", j, result["H"][j]);
        end do;
        for j from 2 to result["Ng"] do
            printf("g[%d] = %a\n", j, result["g"][j]);
        end do;
        for j from 1 to result["Np"] do
            printf("p[%d] = %a\n", j, result["p"][j]);
        end do;
        printf("P(rho) = %a + O(rho^%d)\n",
               result["P_series"], result["Np"]+1);
        printf("Period interpretation requires an independently known center.\n");
        return NULL;
    end proc;
    # 基础分析：显式记录输入、参数条件、证明依据和计算身份。

    # 只记录当前表达式仍可见的分母；输入形成前已经约掉的因子无法恢复。
    VisibleDenominators := proc(expressions::list)
        local result, e, q;
        result := {};
        for e in expressions do
            q := denom(e);
            if not type(q, constant) then result := result union {q}; end if;
            q := denom(normal(e));
            if not type(q, constant) then result := result union {q}; end if;
        end do;
        return [op(result)];
    end proc;

    InputInfo := proc(X, Y, x::name, y::name)
        local parameters, floats, p;
        parameters := remove(type, indets({X,Y}, name) minus {x,y}, constant);
        floats := indets({X,Y}, float);
        return table([
            "X"=X, "Y"=Y, "x"=x, "y"=y,
            "input_degree"=max(degree(X,{x,y}),degree(Y,{x,y})),
            "parameters"=[op(parameters)],
            "active_assumptions"=[op({seq(op(getassumptions(p)),p in parameters)})],
            "exactness"=`if`(floats={}, "exact", "floating"),
            "floating_values"=[op(floats)], "Digits"=Digits,
            "denominator_nonzero"=VisibleDenominators([X,Y]),
            "denominator_scope"="visible expressions only; cancelled input factors cannot be recovered",
            "real_parameter_assumption"=[seq(p::real,p in parameters)],
            "program_version"=Version(), "kernel_version"=kernelopts(version),
            "polar_convention"="x=r*cos(theta), y=r*sin(theta), H0=-1",
            "focus_convention"="g[n]=u[n](2*Pi)/(2*Pi); reverse-angle return",
            "period_convention"="-int(1/H,theta=0..2*Pi); true period only at a center"
        ]);
    end proc;

    RelationPolynomial := proc(e)
        if type(e, equation) then return lhs(e)-rhs(e); end if;
        if type(e, algebraic) then return e; end if;
        error "A relation must be an equation or an expression equal to zero";
    end proc;

    ReduceCondition := proc(f, basis::list, monomial_order)
        local q, numerator, denominator;
        q := normal(f);
        numerator := numer(q);
        denominator := denom(q);
        if basis<>[] then
            numerator := Groebner:-NormalForm(numerator,basis,monomial_order);
            denominator := Groebner:-NormalForm(denominator,basis,monomial_order);
        end if;
        return [numerator,denominator];
    end proc;

    # zero: 有恒等式/理想成员依据；nonzero: 常数或明确非零条件有依据。
    # 非零余式通常仍会在某些参数点消失，必须保留为 unknown。
    # 这里只做理想约化，不声称求出了根理想、全部实分支或非零条件饱和。
    ZeroStatus := proc(f, {relations::list:=[], parameters::list:=[],
                           nonzero::list:=[]})
        local rels, relation_names, parameter_set, basis, monomial_order,
              e, denoms, pair, reduced, result, q, proven, approximate;
        parameter_set := {op(parameters)};
        if nops(parameter_set)<>nops(parameters)
           or not andmap(type,parameters,name) then
            error "parameters must be a list of distinct names";
        end if;
        rels := [seq(RelationPolynomial(e), e in relations)];
        relation_names := remove(type,indets(rels,name),constant);
        if not relation_names subset parameter_set then
            error "parameters must include every symbolic name in the relations";
        end if;
        denoms := VisibleDenominators([f,op(rels),op(nonzero)]);
        approximate := indets([f,op(rels),op(nonzero)],float)<>{};
        basis := [];
        monomial_order := tdeg(op(parameters));
        if rels<>[] then
            if parameters=[] then
                if ormap(e->not evalb(normal(e)=0),rels) then basis := [1]; end if;
            else
                if not andmap(e->type(numer(normal(e)),polynom(anything,parameter_set))
                              and type(denom(normal(e)),polynom(anything,parameter_set)),
                              [f,op(rels),op(nonzero)]) then
                    error "Ideal reduction requires rational polynomial expressions in parameters";
                end if;
                basis := Groebner:-Basis([seq(numer(normal(e)),e in rels)],monomial_order);
            end if;
        end if;
        result := table([
            "status"="unknown", "original"=f, "reduced"=f,
            "relations"=relations, "parameters"=parameters,
            "nonzero"=nonzero, "denominator_nonzero"=denoms,
            "basis"=basis, "domain_status"="not_certified_nonempty",
            "coefficient_domain"="characteristic zero; ideal implication valid over real or complex parameters",
            "evidence"="A nonzero remainder is not a pointwise nonzero certificate",
            "limitations"="No radical membership, branch solving, or saturation is claimed"
        ]);
        if basis=[1] then
            result["domain_status"] := "inconsistent";
            result["evidence"] := "Relations generate the unit ideal; no admissible parameter point";
            return result;
        end if;
        for e in [op(denoms),op(nonzero)] do
            pair := ReduceCondition(e,basis,monomial_order);
            if evalb(pair[1]=0) or evalb(pair[2]=0) then
                result["domain_status"] := "inconsistent";
                result["evidence"] := "A required nonzero expression vanishes modulo the relations";
                return result;
            end if;
        end do;
        pair := ReduceCondition(f,basis,monomial_order);
        if evalb(pair[2]=0) then
            result["domain_status"] := "inconsistent";
            result["evidence"] := "The expression denominator vanishes on the constrained domain";
            return result;
        end if;
        reduced := normal(pair[1]/pair[2]);
        result["reduced"] := reduced;
        if approximate then
            result["evidence"] := "Floating input: no exact algebraic zero/nonzero certificate is issued";
            return result;
        end if;
        if evalb(pair[1]=0) then
            result["status"] := "zero";
            result["evidence"] := `if`(basis=[],"Exact expression identity","Numerator belongs to the supplied polynomial ideal, on the recorded denominator domain");
            return result;
        end if;
        proven := false;
        if type(reduced,constant) then proven := evalb(is(reduced=0)=false); end if;
        for e in nonzero do
            q := ReduceCondition(e,basis,monomial_order);
            q := normal(q[1]/q[2]);
            if evalb(normal(reduced-q)=0) or evalb(normal(reduced+q)=0) then
                proven := true;
            end if;
        end do;
        if proven then
            result["status"] := "nonzero";
            result["evidence"] := "Exact nonzero constant or an explicitly supplied nonzero condition";
        end if;
        return result;
    end proc;

    AnalyzeConditions := proc(result::table, relations::list, parameters::list,
                             {nonzero::list:=[]})
        local g, p, n, domain_nonzero;
        g := table(); p := table();
        # 即使某个输出已恒等为零，也必须保留原方程的可见分母条件。
        domain_nonzero := nonzero;
        if assigned(result["input"]) then
            domain_nonzero := [op({op(nonzero),op(result["input"]["denominator_nonzero"])})];
        end if;
        for n from 2 to result["Ng"] do
            g[n] := ZeroStatus(result["g"][n],':-relations'=relations,
                               ':-parameters'=parameters,':-nonzero'=domain_nonzero);
        end do;
        for n from 0 to result["Np"] do
            p[n] := ZeroStatus(result["p"][n],':-relations'=relations,
                               ':-parameters'=parameters,':-nonzero'=domain_nonzero);
        end do;
        return table([
            "g"=g, "p"=p, "relations"=relations,
            "parameters"=parameters, "nonzero"=domain_nonzero,
            "raw_result_unchanged"=true,
            "interpretation"="Ideal implications only; a finite zero focus list is not a center proof",
            "domain"="real or complex algebraic parameter constraints; real feasibility not decided"
        ]);
    end proc;

    # 只给出有独立依据的中心充分条件；未知参数按实参数前提明确列出。
    # 即使全部已算 g 为零，也不因此认定中心。
    CenterCertificate := proc(result::table, {first_integral:=0})
        local input, X, Y, x, y, parameters, assumptions, coefficients,
              c, real_check, out, divergence, rx, ry, invariant,
              gradient, h11, determinant, positive, negative;
        if not assigned(result["input"]) then
            error "Result has no input record; recompute with the current module";
        end if;
        input := result["input"];
        X := input["X"]; Y := input["Y"];
        x := input["x"]; y := input["y"];
        parameters := input["parameters"];
        assumptions := [seq(c::real,c in parameters)];
        out := table(["status"="unknown", "method"="none",
                      "evidence"="No supported center certificate was verified",
                      "assumptions"=assumptions, "input"=input,
                      "denominator_nonzero"=input["denominator_nonzero"],
                      "real_coefficients"=false,
                      "finite_focus_quantities_are_not_a_proof"=true]);
        # InputInfo是公开的元数据构造器；中心定理不能仅信任手工表的标准性。
        if x=y or not type(X,polynom(anything,{x,y}))
           or not type(Y,polynom(anything,{x,y})) then
            out["evidence"] := "A standard polynomial vector field is required";
            return out;
        end if;
        if not evalb(Clean(eval(X,{x=0,y=0}))=0)
           or not evalb(Clean(eval(Y,{x=0,y=0}))=0)
           or not evalb(Clean(eval(diff(X,x),{x=0,y=0}))=0)
           or not evalb(Clean(eval(diff(X,y),{x=0,y=0})-1)=0)
           or not evalb(Clean(eval(diff(Y,x),{x=0,y=0})+1)=0)
           or not evalb(Clean(eval(diff(Y,y),{x=0,y=0}))=0) then
            out["evidence"] := "The equilibrium and fixed elliptic linear-part premises were not verified";
            return out;
        end if;
        if input["exactness"]<>"exact" then
            out["evidence"] := "Floating coefficients do not supply an exact center certificate";
            return out;
        end if;
        coefficients := [coeffs(expand(X),{x,y}),coeffs(expand(Y),{x,y})];
        for c in coefficients do
            real_check := is(c,'real') assuming op(assumptions);
            if not evalb(real_check=true) then
                out["evidence"] := "Real coefficients could not be verified under the recorded real-parameter assumptions";
                return out;
            end if;
        end do;
        out["real_coefficients"] := true;
        divergence := Clean(diff(X,x)+diff(Y,y));
        if evalb(divergence=0) then
            out["status"] := "proved";
            out["method"] := "Hamiltonian";
            out["evidence"] := "Zero divergence gives a local analytic Hamiltonian; the fixed linear part gives a positive definite quadratic Hamiltonian";
            out["divergence"] := divergence;
            return out;
        end if;
        rx := Clean(subs(y=-y,X)+X);
        ry := Clean(subs(y=-y,Y)-Y);
        if evalb(rx=0) and evalb(ry=0) then
            out["status"] := "proved";
            out["method"] := "reversible_y_reflection";
            out["evidence"] := "Analytic real vector field is reversible under (x,y)->(x,-y), with the fixed nondegenerate elliptic linear part";
            out["reversibility_residual"] := [rx,ry];
            return out;
        end if;
        if evalb(first_integral=0) then return out; end if;
        if indets(first_integral,float)<>{} then
            out["evidence"] := "The proposed first integral contains floating values";
            return out;
        end if;
        # 当前充分条件核验覆盖多项式第一积分；不把未核验的解析性当作事实。
        if not type(first_integral,polynom(anything,{x,y})) then
            out["evidence"] := "This certificate checks polynomial first integrals; analyticity of a general expression was not certified";
            return out;
        end if;
        real_check := is(first_integral,'real') assuming op(assumptions),x::real,y::real;
        if not evalb(real_check=true) then
            out["evidence"] := "The proposed first integral was not verified to be real under recorded assumptions";
            return out;
        end if;
        invariant := Clean(diff(first_integral,x)*X+diff(first_integral,y)*Y);
        gradient := [Clean(eval(diff(first_integral,x),{x=0,y=0})),
                     Clean(eval(diff(first_integral,y),{x=0,y=0}))];
        h11 := Clean(eval(diff(first_integral,x,x),{x=0,y=0}));
        determinant := Clean(h11*eval(diff(first_integral,y,y),{x=0,y=0})
                             -eval(diff(first_integral,x,y),{x=0,y=0})^2);
        positive := is(h11>0) assuming op(assumptions);
        negative := is(h11<0) assuming op(assumptions);
        real_check := is(determinant>0) assuming op(assumptions);
        out["first_integral"] := first_integral;
        out["first_integral_residual"] := invariant;
        out["gradient_at_origin"] := gradient;
        out["hessian_principal_minors"] := [h11,determinant];
        if evalb(invariant=0) and gradient=[0,0]
           and evalb(real_check=true)
           and (evalb(positive=true) or evalb(negative=true)) then
            out["status"] := "proved";
            out["method"] := "first_integral_extremum";
            out["evidence"] := "Verified real first integral with zero gradient and definite Hessian at the origin";
        else
            out["evidence"] := "The first-integral identity and a nondegenerate real extremum were not both verified";
        end if;
        return out;
    end proc;

    # source_file 是用户显式指定的文件，不能据此声称它就是已经加载的源码。
    # Maple 2025 内置 FileTools:-Hash 支持 SHA1，而不支持 SHA256。
    # 此指纹用于辨识文件版本，不用于安全签名；发行附件可在外部另算 SHA256。
    RunInfo := proc(result::table, {source_file::string:="", revision::string:=""})
        local info;
        info := table([
            "program_version"="not_recorded", "kernel_version"="not_recorded",
            "Digits"="not_recorded", "revision"=revision,
            "inspection_kernel_version"=kernelopts(version),
            "orders"=["Ng"=result["Ng"],"Np"=result["Np"],
                      "u_order"=result["u_order"],"A_order"=result["A_order"]],
            "direction"="theta increases from 0 to 2*Pi; focus map is reverse-angle; period density carries a minus sign",
            "source_file"=source_file, "source_digest"="",
            "source_digest_algorithm"="", "source_identity_status"="not_recorded",
            "source_identity_scope"="A user-specified file digest does not prove which source was loaded"
        ]);
        if assigned(result["input"]) then
            info["input"] := result["input"];
            info["program_version"] := result["input"]["program_version"];
            info["kernel_version"] := result["input"]["kernel_version"];
            info["Digits"] := result["input"]["Digits"];
        end if;
        info["method"] := result["method"];
        info["optimize"] := result["optimize"];
        info["factor_orbits"] := result["factor_orbits"];
        info["compact"] := result["compact"];
        info["orbit_factor"] := result["orbit_factor"];
        if assigned(result["jet"]) then info["jet"] := result["jet"]; end if;
        if source_file<>"" then
            if not FileTools:-Exists(source_file) or FileTools:-IsDirectory(source_file) then
                error "source_file must name an existing local file";
            end if;
            info["source_digest"] := FileTools:-Hash(source_file,method="sha1");
            info["source_digest_algorithm"] := "sha1";
            info["source_identity_status"] := "specified_file_hashed";
        end if;
        return info;
    end proc;

    # 参数局部展开贯穿 Compute 的 jet 内核；不是算完全参数通式后再截断。
    # epsilon 是本过程私有符号；返回后通过 ans["epsilon"] 查询。
    ParameterJet := proc(X, Y, x::name, y::name, Ng::nonnegint, Np::nonnegint,
                         base::set, direction::set, K::nonnegint:=1)
        local epsilon, base_table, direction_table, substitutions, e, p,
              result, g_coefficients, p_coefficients, n, k;
        base_table := table(); direction_table := table();
        for e in base do
            if not type(e,equation) or not type(lhs(e),name)
               or member(lhs(e),{x,y}) then
                error "base must contain parameter=value equations and no state variables";
            end if;
            if assigned(base_table[lhs(e)]) then error "Duplicate base parameter"; end if;
            base_table[lhs(e)] := rhs(e);
        end do;
        for e in direction do
            if not type(e,equation) or not type(lhs(e),name)
               or member(lhs(e),{x,y}) then
                error "direction must contain parameter=value equations and no state variables";
            end if;
            if assigned(direction_table[lhs(e)]) then error "Duplicate direction parameter"; end if;
            if not assigned(base_table[lhs(e)]) then
                error "Every direction parameter must have an explicit base value";
            end if;
            direction_table[lhs(e)] := rhs(e);
        end do;
        substitutions := {};
        for p in [indices(base_table,nolist)] do
            if has(base_table[p],{x,y}) then error "Base values must not depend on state variables"; end if;
            if assigned(direction_table[p]) then
                if has(direction_table[p],{x,y}) then error "Directions must not depend on state variables"; end if;
                substitutions := substitutions union {p=base_table[p]+epsilon*direction_table[p]};
            else
                substitutions := substitutions union {p=base_table[p]};
            end if;
        end do;
        result := Compute(subs(substitutions,X),subs(substitutions,Y),x,y,Ng,Np,
                          jet=[epsilon,K]);
        g_coefficients := table(); p_coefficients := table();
        for n from 2 to Ng do
            for k from 0 to K do
                g_coefficients[n,k] := coeff(expand(result["g"][n]),epsilon,k);
            end do;
        end do;
        for n from 0 to Np do
            for k from 0 to K do
                p_coefficients[n,k] := coeff(expand(result["p"][n]),epsilon,k);
            end do;
        end do;
        result["epsilon"] := epsilon;
        result["parameter_jet"] := table([
            "base"=base,"direction"=direction,"order"=K,
            "original_input"=InputInfo(X,Y,x,y),"substitutions"=substitutions,
            "epsilon"=epsilon,"g_coefficients"=g_coefficients,
            "p_coefficients"=p_coefficients,
            "meaning"="Taylor coefficients along the specified parameter direction; not a full-parameter formula",
            "center_boundary"="Parameter directions need not stay in a center family; period coefficients may describe formal time only"
        ]);
        return result;
    end proc;

    # 仅按x,y总次数取齐次分量；参数不计入次数。
    HomogeneousPart := proc(f, x::name, y::name, m::nonnegint)
        local e, j;
        e := expand(f);
        return add(coeff(coeff(e,x,m-j),y,j)*x^(m-j)*y^j,j=0..m);
    end proc;

    # 独立笛卡尔内核：N是第一积分/残差的最高总次数，不是g的下标。
    # V=x^2+y^2+sum(V_m)，X*V_x+Y*V_y=sum(L[k]*(x^2+y^2)^(k+1))。
    # 偶数m固定V_m中x^m系数为0，消除径向核；L不是用户定义的完整g。
    Lyapunov := proc(X, Y, x::name, y::name, N::posint)
        local XX, YY, d, Fx, Fy, Vc, L, V, target, residual,
              m, j, i, k, W, M, rhs, sol, even_degree, size_m;
        if x=y or N<2 then
            error "Distinct state names and total degree N>=2 are required";
        end if;
        XX := expand(X); YY := expand(Y);
        if not type(XX,polynom(anything,{x,y}))
           or not type(YY,polynom(anything,{x,y})) then
            error "X and Y must be polynomials in the state variables";
        end if;
        if not evalb(simplify(eval(XX,{x=0,y=0}))=0)
           or not evalb(simplify(eval(YY,{x=0,y=0}))=0)
           or not evalb(simplify(eval(diff(XX,x),{x=0,y=0}))=0)
           or not evalb(simplify(eval(diff(XX,y),{x=0,y=0})-1)=0)
           or not evalb(simplify(eval(diff(YY,x),{x=0,y=0})+1)=0)
           or not evalb(simplify(eval(diff(YY,y),{x=0,y=0}))=0) then
            error "The equilibrium and linear matrix must be (0,0), [[0,1],[-1,0]]";
        end if;
        d := max(degree(XX,{x,y}),degree(YY,{x,y}));
        Fx := Array(1..max(1,d),fill=0);
        Fy := Array(1..max(1,d),fill=0);
        for j from 2 to d do
            Fx[j] := HomogeneousPart(XX,x,y,j);
            Fy[j] := HomogeneousPart(YY,x,y,j);
        end do;
        Vc := Array(2..N,fill=0); Vc[2] := x^2+y^2;
        L := table(); V := Vc[2]; target := 0;
        for m from 3 to N do
            W := expand(add(diff(Vc[m-j+1],x)*Fx[j]
                      +diff(Vc[m-j+1],y)*Fy[j],j=2..min(d,m-1)));
            even_degree := irem(m,2)=0;
            size_m := m+1;
            if even_degree then size_m := m+2; end if;
            M := Matrix(size_m,size_m,fill=0);
            rhs := Vector(size_m,fill=0);
            for i from 0 to m do
                rhs[i+1] := -coeff(coeff(W,x,m-i),y,i);
                if i<m then M[i+2,i+1] := m-i; end if;
                if i>0 then M[i,i+1] := -i; end if;
                if even_degree and irem(i,2)=0 then
                    M[i+1,size_m] := -binomial(m/2,i/2);
                end if;
            end do;
            if even_degree then M[size_m,1] := 1; end if;
            sol := LinearAlgebra:-LinearSolve(M,rhs);
            Vc[m] := expand(add(sol[i+1]*x^(m-i)*y^i,i=0..m));
            V := V+Vc[m];
            if even_degree then
                k := m/2-1; L[k] := simplify(sol[size_m]);
                target := target+L[k]*(x^2+y^2)^(k+1);
            end if;
        end do;
        residual := expand(diff(V,x)*XX+diff(V,y)*YY-target);
        residual := simplify(add(HomogeneousPart(residual,x,y,m),m=0..N));
        return table(["degree"=N,"V"=expand(V),"V_components"=Vc,
            "L"=L,"target_series"=expand(target),"residual"=residual,
            "input"=InputInfo(XX,YY,x,y),
            "normalization"="V2=x^2+y^2; coefficient of x^m in Vm is zero for even m",
            "interpretation"="Independent Cartesian Lyapunov quantities; not the complete g[n]",
            "first_nonzero_relation"="When L1..L(k-1)=0, g[2*k+1]=-L[k]/2 in the prescribed angle direction"]);
    end proc;

    # 按原始非线性系数的权(d-1)核验坐标缩放律。
    # 这是显式运行的交叉检查，不在每次Compute中重复计算。
    ScalingCheck := proc(result::table, scale::algebraic)
        local input, x, y, Xs, Ys, scaled, gres, pres, n, m, status;
        input := result["input"]; x := input["x"]; y := input["y"];
        if evalb(scale=0) or has(scale,{x,y,result["theta"],result["r"],result["rho"]}) then
            error "Use a nonzero scaling factor independent of the state and polar variables";
        end if;
        Xs := add(HomogeneousPart(input["X"],x,y,m)*scale^(m-1),m=1..result["degree"]);
        Ys := add(HomogeneousPart(input["Y"],x,y,m)*scale^(m-1),m=1..result["degree"]);
        scaled := Compute(Xs,Ys,x,y,result["Ng"],result["Np"],
                          method=result["method"],jet=result["jet"]);
        gres := table(); pres := table(); status := "verified";
        for n from 2 to result["Ng"] do
            gres[n] := Clean(scaled["g"][n]-scale^(n-1)*result["g"][n],result["jet"]);
            if not evalb(gres[n]=0) then status := "undetermined"; end if;
        end do;
        for n from 0 to result["Np"] do
            pres[n] := Clean(scaled["p"][n]-scale^n*result["p"][n],result["jet"]);
            if not evalb(pres[n]=0) then status := "undetermined"; end if;
        end do;
        if input["exactness"]<>"exact" or hastype(scale,float) then
            status := "approximate_check";
        end if;
        return table(["scale"=scale,"scale_nonzero"=(scale<>0),
            "status"=status,"g_residual"=gres,"p_residual"=pres,
            "g_weight"="n-1","p_weight"="n",
            "interpretation"="coordinate scaling of independent polynomial coefficients; does not assign one weight to a reused mixed-degree parameter"]);
    end proc;

    # 返回映射与周期研究：所有新量单独保存，不覆盖原始g/p。

    # 显式传入展开变量，兼容来自另一模块加载实例的结果私有符号。
    MathCut := proc(f, variable::name, N::nonnegint, jet::list := [])
        local e, j;
        e := expand(JetCut(f,jet));
        return JetCut(add(coeff(e,variable,j)*variable^j,j=0..N),jet);
    end proc;

    # Horner 复合，每一步截断；不先构造完整高次复合式。
    MathCompose := proc(f, g, variable::name, N::nonnegint, jet::list := [])
        local out, j, e;
        out := 0;
        e := expand(JetCut(f,jet));
        for j from N by -1 to 0 do
            out := MathCut(out*g,variable,N,jet)+coeff(e,variable,j);
        end do;
        return MathCut(out,variable,N,jet);
    end proc;

    # 只用于导数为 1 的近恒等形式级数，不涉及数值求逆或分支选择。
    MathInverse := proc(f, variable::name, N::posint, jet::list := [])
        local inverse, n;
        if not evalb(Clean(coeff(expand(f),variable,0))=0)
           or not evalb(Clean(coeff(expand(f),variable,1))=1) then
            error "Formal inversion requires constant coefficient 0 and linear coefficient 1";
        end if;
        inverse := variable;
        for n from 2 to N do
            inverse := inverse-Clean(coeff(MathCompose(f,inverse,variable,n,jet),variable,n),jet)*variable^n;
        end do;
        return MathCut(inverse,variable,N,jet);
    end proc;

    # 只为无参数的精确实常数判定符号，不借用工作表对参数的隐藏 assume 状态。
    MathConstantSign := proc(f)
        local c;
        c := Clean(f);
        if evalb(c=0) then return "zero"; end if;
        if hastype(c,float) then return "unknown"; end if;
        if type(c,realcons) then
            if is(c>0)=true then return "positive"; end if;
            if is(c<0)=true then return "negative"; end if;
        end if;
        return "unknown";
    end proc;

    # 半周图 S 与带符号半周图 J。J o J 是同一原截面的全周返回图。
    # 不约化、不覆盖用户的完整焦点量 g[n]。
    HalfReturn := proc(result::table)
        local N, angle, radius, S, J, full, composed, residual,
              dependencies, n, status, jet;
        N := result["u_order"];
        jet := [];
        if assigned(result["jet"]) then jet := result["jet"]; end if;
        angle := result["theta"];
        radius := result["rho"];
        S := Clean(add(eval(result["u"][n],angle=Pi)*radius^n,n=1..N),jet);
        J := -S;
        full := Clean(add(eval(result["u"][n],angle=2*Pi)*radius^n,n=1..N),jet);
        composed := Clean(MathCompose(J,J,radius,N,jet),jet);
        residual := Clean(composed-full,jet);
        status := "undetermined";
        if evalb(residual=0) then status := "verified_to_order"; end if;
        dependencies := table();
        if N>=2 then
            dependencies["g2_residual"] := Clean(coeff(full-radius,radius,2)/(2*Pi));
        end if;
        if N>=4 then
            dependencies["s2"] := Clean(coeff(S,radius,2));
            dependencies["g4_minus_s2_g3_over_2"] := Clean(
                coeff(full-radius,radius,4)/(2*Pi)
                -coeff(S,radius,2)*coeff(full-radius,radius,3)/(4*Pi),jet);
        end if;
        return table([
            "rho"=radius,"order"=N,"jet"=jet,"S_series"=S,"J_series"=J,
            "full_return_series"=full,"half_composed_series"=composed,
            "composition_residual"=residual,"composition_status"=status,
            "dependencies"=dependencies,
            "interpretation"="signed half-return identity; finite-order agreement is not a center proof"
        ]);
    end proc;

    # 当前 g 对应逆时间角向图；正时间图是它的形式逆，而非逐阶取负。
    ForwardReturn := proc(result::table)
        local half, N, radius, inverse, residual, g, n, c, sign_status,
              stability, leading_order, candidate_order, candidate_coefficient, jet;
        half := HalfReturn(result);
        N := half["order"];
        radius := half["rho"];
        jet := half["jet"];
        inverse := Clean(MathInverse(half["full_return_series"],radius,N,jet),jet);
        residual := Clean(MathCompose(half["full_return_series"],inverse,radius,N,jet)-radius,jet);
        g := table();
        for n from 2 to N do
            g[n] := Clean(coeff(inverse-radius,radius,n)/(2*Pi),jet);
        end do;
        stability := "unknown";
        leading_order := "undetermined";
        candidate_order := "none_to_computed_order";
        candidate_coefficient := 0;
        for n from 2 to N do
            c := Clean(coeff(inverse-radius,radius,n));
            sign_status := MathConstantSign(c);
            if sign_status="zero" then next; end if;
            candidate_order := n;
            candidate_coefficient := c;
            if sign_status="negative" then
                stability := "stable";
                leading_order := n;
            elif sign_status="positive" then
                stability := "unstable";
                leading_order := n;
            end if;
            break;
        end do;
        return table([
            "rho"=radius,"order"=N,"jet"=jet,"return_series"=inverse,"g"=g,
            "inverse_residual"=residual,"stability"=stability,
            "first_nonzero_order"=leading_order,
            "leading_candidate_order"=candidate_order,
            "leading_candidate_coefficient"=candidate_coefficient,
            "direction"="positive physical time; theta decreases by 2*Pi",
            "interpretation"="local formal inverse; all computed zero coefficients do not prove a center"
        ]);
    end proc;

    # 中心的自然振幅 alpha=(rho-J(rho))/2，将半周图线性化为 alpha -> -alpha。
    # 默认必须有独立中心证据；formal=true 只做代数转换，不把有限消失误称中心。
    CenterPeriod := proc(result::table, N::integer := -1,
                         {first_integral := 0, formal::truefalse := false})
        local requested, half, radius, certificate, transform, inverse,
              P, P_alpha, tau, odd_residual, even_series, inverse_residual,
              n, k, status, interpretation, jet;
        requested := N;
        if requested=-1 then requested := result["Np"]; end if;
        if requested<0 or requested>result["Np"] then
            error "Requested period order must lie between 0 and the computed Np";
        end if;
        if result["u_order"]<max(1,requested) then
            error "Insufficient radial order for amplitude conversion; recompute with Ng at least %1", max(2,requested);
        end if;
        # 全局限定的关键字标签避免与本过程同名关键字参数混淆。
        certificate := CenterCertificate(result,':-first_integral'=first_integral);
        status := certificate["status"];
        if status<>"proved" and not formal then
            error "A verified center certificate is required; formal=true only permits an explicitly unproved transformation";
        end if;
        half := HalfReturn(result);
        jet := half["jet"];
        radius := result["rho"];
        transform := MathCut((radius-half["J_series"])/2,radius,max(1,requested),jet);
        inverse := MathInverse(transform,radius,max(1,requested),jet);
        inverse_residual := Clean(MathCompose(transform,inverse,radius,max(1,requested),jet)-radius,jet);
        P := add(result["p"][n]*radius^n,n=0..requested);
        P_alpha := Clean(subs(radius=alpha,MathCompose(P,inverse,radius,requested,jet)),jet);
        tau := table();
        for k from 1 to iquo(requested,2) do
            tau[k] := Clean(coeff(P_alpha,alpha,2*k));
        end do;
        odd_residual := Clean(add(coeff(P_alpha,alpha,2*k+1)*alpha^(2*k+1),
                                  k=0..iquo(requested-1,2)));
        even_series := Clean(coeff(P_alpha,alpha,0)
                         +add(tau[k]*alpha^(2*k),k=1..iquo(requested,2)));
        interpretation := "physical center period";
        if status<>"proved" then
            status := "formal_only";
            interpretation := "unproved amplitude transformation of formal time; not a center certificate";
        end if;
        return table([
            "rho"=radius,"alpha"=alpha,"order"=requested,"jet"=jet,
            "alpha_series"=transform,"rho_of_alpha"=subs(radius=alpha,inverse),
            "inverse_residual"=inverse_residual,"P_alpha_series"=P_alpha,
            "P_even_series"=even_series,"tau"=tau,"tau_order"=iquo(requested,2),
            "odd_residual"=odd_residual,"center_status"=status,
            "center_certificate"=certificate,"half_return"=half,
            "interpretation"=interpretation
        ]);
    end proc;

    # 截断临界周期候选：P'(alpha)/(2*alpha)=sum(k*tau[k]*z^(k-1)), z=alpha^2>0。
    # 正根仅是截断候选；真实临界周期需要余项控制及分岔定理。
    CriticalPeriods := proc(period_data::table, parameters::list := [])
        local K, polynomial, roots, candidates, root, sign_status, positivity,
              jacobian, generic_rank, i, j, k, degree_status, scope,
              leading_coefficient, branch_conditions, jet, rank_scope,
              interpretation;
        if not assigned(period_data["tau_order"]) or not assigned(period_data["tau"]) then
            error "Use the table returned by CenterPeriod";
        end if;
        if not andmap(type,parameters,name) or nops({op(parameters)})<>nops(parameters) then
            error "Parameters must be a list of distinct unassigned names";
        end if;
        K := period_data["tau_order"];
        jet := [];
        if assigned(period_data["jet"]) then jet := period_data["jet"]; end if;
        if not assigned(period_data["odd_residual"]) then
            error "Use CenterPeriod output, including its odd_residual";
        end if;
        if not evalb(Clean(period_data["odd_residual"],jet)=0) then
            error "The amplitude series has odd terms; a polynomial in z=alpha^2 does not describe its full derivative";
        end if;
        polynomial := Clean(add(k*period_data["tau"][k]*z^(k-1),k=1..K));
        roots := [];
        leading_coefficient := 0;
        branch_conditions := [];
        degree_status := "computed_truncation";
        scope := "generic parameter formula; degenerate leading coefficients and denominators require separate branch analysis";
        if evalb(polynomial=0) then
            degree_status := "no_information_beyond_computed_order";
        elif degree(polynomial,z)<=4 then
            leading_coefficient := coeff(polynomial,z,degree(polynomial,z));
            if MathConstantSign(leading_coefficient)="unknown" then
                branch_conditions := [leading_coefficient<>0];
            end if;
            roots := [solve(polynomial=0,z)];
        else
            leading_coefficient := coeff(polynomial,z,degree(polynomial,z));
            degree_status := "equation_returned_without_explicit_roots";
        end if;
        candidates := [];
        for root in roots do
            sign_status := MathConstantSign(root);
            positivity := "undetermined";
            if sign_status="positive" then positivity := "positive";
            elif member(sign_status,{"zero","negative"}) then positivity := "nonpositive";
            elif type(root,complexcons) and not type(root,realcons) then positivity := "not_real";
            end if;
            candidates := [op(candidates),table(["z"=root,"positivity"=positivity])];
        end do;
        jacobian := Matrix(K,nops(parameters),
                           (i,j)->diff(period_data["tau"][i],parameters[j]));
        generic_rank := LinearAlgebra:-Rank(jacobian);
        rank_scope := "generic algebraic rank; may decrease on special parameter branches";
        interpretation := "truncated candidates only; parameter directions must preserve the independently verified center family";
        if jet<>[] then
            rank_scope := "rank of the truncated parameter expression only; not a certified rank of the untruncated family";
            interpretation := "candidates from a truncated radial and parameter expansion only; omitted parameter powers are unknown, and directions must preserve a verified center family";
        end if;
        return table([
            "z"=z,"jet"=jet,"critical_polynomial"=polynomial,
            "candidate_equation"=(polynomial=0),"candidate_domain"=(z>0),
            "candidates"=candidates,"degree_status"=degree_status,
            "leading_coefficient"=leading_coefficient,"generic_branch_conditions"=branch_conditions,
            "parameters"=parameters,"jacobian"=jacobian,"generic_rank"=generic_rank,
            "rank_scope"=rank_scope,
            "solution_scope"=scope,"center_status"=period_data["center_status"],
            "interpretation"=interpretation
        ]);
    end proc;

end module:
