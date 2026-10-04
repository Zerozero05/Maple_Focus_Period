# ======================================================================
# 标准平面多项式系统：焦点量与周期系数
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
    export Compute, PolarData, Show, Version;
    local Clean, Cut, IntegralFromZero, Reciprocal, RadialSolution,
          PeriodSeries, theta, r, rho;

    # theta、r、rho 是模块私有符号，避免与输入中的同名参数混淆。
    # 需要代入时，请使用 ans["theta"]、ans["r"]、ans["rho"]。

    # 每个发布版本同步维护此版本号、根目录 VERSION 和 CHANGELOG.md。
    Version := proc()
        return "1.0.1";
    end proc;

    Clean := proc(f)
        return simplify(expand(f), trig);
    end proc;

    # 对 rho 多项式精确截断，避免直接展开高次幂产生不需要的项。
    Cut := proc(f, N::nonnegint)
        local e, k;
        e := expand(f);
        return add(coeff(e, rho, k)*rho^k, k=0..N);
    end proc;

    # 从 0 积分到 upper；减去积分常数，保证 u[n](0)=0。
    # 若 Maple 留下未计算的积分，明确停止，不将其冒充为已求出的系数。
    IntegralFromZero := proc(f, upper)
        local primitive;
        if evalb(f=0) then
            return 0;
        end if;
        primitive := int(f, theta);
        if has(primitive, {'int', 'Int'}) then
            error "Maple could not evaluate an integral; integrand: %1", f;
        end if;
        return Clean(eval(primitive, theta=upper)
                     - eval(primitive, theta=0));
    end proc;

    # 第1步：检查原始系统，自动极坐标化并按 r 的幂提取 G、H。
    PolarData := proc(X, Y, x::name, y::name)
        local XX, YY, d, xp, yp, radial, angular, G, H, j;
        if evalb(x=y) then
            error "The two state variables must be distinct";
        end if;
        if not type(X, polynom(anything, {x,y}))
           or not type(Y, polynom(anything, {x,y})) then
            error "X and Y must be polynomials in the state variables";
        end if;
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

        d := max(degree(XX, {x,y}), degree(YY, {x,y}));
        xp := subs({x=r*cos(theta), y=r*sin(theta)}, XX);
        yp := subs({x=r*cos(theta), y=r*sin(theta)}, YY);
        radial := expand(xp*cos(theta) + yp*sin(theta));
        angular := expand((yp*cos(theta) - xp*sin(theta))/r);

        G := Array(0..d, fill=0);
        H := Array(0..d-1, fill=0);
        for j from 0 to d do
            G[j] := Clean(coeff(radial, r, j));
        end do;
        for j from 0 to d-1 do
            H[j] := Clean(coeff(angular, r, j));
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
    Reciprocal := proc(H::Array, d::posint, N::nonnegint)
        local A, n, j;
        A := Array(0..N, fill=0);
        A[0] := 1/H[0];
        for n from 1 to N do
            A[n] := Clean(-add(H[j]*A[n-j], j=1..min(n,d-1))/H[0]);
        end do;
        return A;
    end proc;

    # 第3步：dr/dtheta=sum(R[j]*r^j,j>=2)，
    # R[j]=sum(G[k]*A[j-k],k=2..min(j,d))。
    # r(theta,rho)=rho+sum(u[n](theta)*rho^n,n>=2)。
    # 求 u[n]' 时只需要 u[1]..u[n-1]，因为径向方程从 r^2 开始。
    RadialSolution := proc(G::Array, A::Array, d::posint, N::posint)
        local R, u, du, solution, power, rhs, n, j, k;
        R := Array(0..N, fill=0);
        u := Array(1..N, fill=0);
        du := Array(1..N, fill=0);
        u[1] := 1;
        solution := rho;
        for j from 2 to N do
            R[j] := Clean(add(G[k]*A[j-k], k=2..min(j,d)));
        end do;
        for n from 2 to N do
            power := solution;
            rhs := 0;
            for j from 2 to n do
                power := Cut(power*solution, n);
                rhs := rhs + R[j]*coeff(power, rho, n);
            end do;
            du[n] := Clean(rhs);
            u[n] := IntegralFromZero(du[n], theta);
            solution := solution + u[n]*rho^n;
        end do;
        return table(["R"=R, "u"=u, "du"=du, "r_series"=solution]);
    end proc;

    # 第4步：先将 r(theta,rho) 代入 sum(A[j]*r^j)，再提取 rho 系数。
    # 注意不能直接对 A[n] 积分；A[1]*u[n] 等混合项同样贡献 p[n]。
    # 负号只作用在时间密度上，保持原极坐标与 u[n] 完全不变。
    PeriodSeries := proc(A::Array, solution, N::nonnegint)
        local p, B, rp, power, density, j, n;
        p := Array(0..N, fill=0);
        B := Array(0..N, fill=0);
        rp := Cut(solution, N);
        power := 1;
        density := A[0];
        for j from 1 to N do
            power := Cut(power*rp, N);
            density := density + A[j]*power;
        end do;
        density := -Cut(density, N);
        for n from 0 to N do
            B[n] := Clean(coeff(density, rho, n));
            p[n] := IntegralFromZero(B[n], 2*Pi);
        end do;
        return table([
            "p"=p, "B"=B,
            "time_density_series"=add(B[n]*rho^n, n=0..N),
            "P_series"=add(p[n]*rho^n, n=0..N)
        ]);
    end proc;

    # 主接口：两个最高阶独立指定；这里只计算实际需要的中间阶数。
    # u 需要到 max(Ng,Np,1)，A 需要到 max(Np,N_u-2,0)。
    Compute := proc(X, Y, x::name, y::name,
                    Ng::nonnegint, Np::nonnegint)
        local result, d, Nu, NA, A, radial, period, g, n;
        if Ng=1 then
            error "Use Ng=0 to omit focus quantities, or Ng>=2";
        end if;
        result := PolarData(X, Y, x, y);
        d := result["degree"];
        Nu := max(1, Ng, Np);
        NA := max(0, Np, Nu-2);
        A := Reciprocal(result["H"], d, NA);
        radial := RadialSolution(result["G"], A, d, Nu);
        period := PeriodSeries(A, radial["r_series"], Np);
        g := table();
        for n from 2 to Ng do
            g[n] := Clean(eval(radial["u"][n], theta=2*Pi)/(2*Pi));
        end do;

        result["Ng"] := Ng;
        result["Np"] := Np;
        result["u_order"] := Nu;
        result["A_order"] := NA;
        result["A"] := A;
        result["R"] := radial["R"];
        result["u"] := radial["u"];
        result["du"] := radial["du"];
        result["g"] := g;
        result["p"] := period["p"];
        result["B"] := period["B"];
        result["r_series"] := radial["r_series"];
        result["P_series"] := period["P_series"];
        result["time_density_series"] := period["time_density_series"];
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
end module:
