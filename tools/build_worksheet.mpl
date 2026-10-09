# ======================================================================
# 生成原生 Maple .mw 计算工作表（不是给 .mpl 更换扩展名）。
# 在仓库根目录执行：read "tools/build_worksheet.mpl":
# 产物：examples/FocusPeriod_Worksheet.mw。
# 使用官方 DocumentTools:-Layout 和 Worksheet:-ToString API。
# ContentToString 需要 Standard Worksheet GUI；ToString 可在命令行使用。
# Windows 路径字符串使用 / 或 \\；单反斜杠会被当作转义符处理。
# ======================================================================

BuildWorksheet := proc(outputFile::string,
                       repoRoot::string := "D:/Maple/Maple_Focus_Period")
    local L, Code, Text, Heading, worksheet, mwtext, styles;
    L := DocumentTools:-Layout;
    Code := proc(source::string)
        return DocumentTools:-Layout:-Group(
            DocumentTools:-Layout:-Input(
                DocumentTools:-Layout:-Textfield(source,
                    style=MapleInput, layout=Normal,
                    prompt="> ", alignment=left)));
    end proc;
    Text := proc(source::string)
        return DocumentTools:-Layout:-Group(
            DocumentTools:-Layout:-Input(
                DocumentTools:-Layout:-Textfield(source,
                    style=:-Text, layout=Normal, alignment=left)));
    end proc;
    Heading := proc(source::string)
        return DocumentTools:-Layout:-Group(
            DocumentTools:-Layout:-Input(
                DocumentTools:-Layout:-Textfield(source,
                    style=Heading2, layout=Heading2, alignment=left)));
    end proc;

    worksheet := L:-Worksheet(
        L:-Group(L:-Input(L:-Textfield(
            "Maple Focus Period：焦点量与周期系数计算工作表",
            style=Heading1, layout=Heading1, alignment=left))),
        Text("先修改下方 repoRoot 为下载仓库的根目录，再由上至下执行各组输入。此示例执行 restart，建议在独立工作表中运行。核心算法只保存在 src/focus_period.mpl；本工作表用于设置系统、阶数并查看计算结果。"),
        Text("Windows 路径必须正确书写：推荐用正斜杠 /，例如 D:/Maple/Maple_Focus_Period。若使用反斜杠，每个分隔符须写成双反斜杠。不要把资源管理器的单反斜杠路径直接粘入字符串，否则目录分隔符会被按转义字符处理，导致 currentdir、read 失败，随后报 FocusPeriod 不是模块。项目根目录应直接包含 src、examples 等文件夹。"),
        Text("红色的一维 Maple Input 是可执行代码，可以直接运行，无须转换成二维数学公式；普通文字区仅用于说明，不执行。请先执行第 1 块并确认版本显示为 1.0.2，再运行其余输入块；也可点击执行工作表。"),
        Text("固定 x=r*cos(theta)、y=r*sin(theta)，所以 H[0]=-1。g[n] 采用 u[n](2*Pi)/(2*Pi)，n>=2。P 的时间密度取负号，保证 P(0)=2*Pi；只有已独立证明为中心时，P 才是真实周期函数。"),

        Heading("1. 加载核心程序"),
        Code(cat("restart:\n# Windows 路径用 /；若用反斜杠，每个分隔符要写两次。\n# 把下面的示例目录改成你实际的项目根目录（直接包含 src）。\nrepoRoot := ",sprintf("%a",repoRoot),":\n",
                 "currentdir(repoRoot):\nread \"src/focus_period.mpl\":\n",
                 "FocusPeriod:-Version();")),

        Heading("2. 输入原始多项式系统与两个最高阶"),
        Text("以下系统为解析中心，含二次、三次项。用自己的原始 X、Y 替换这两行即可；线性矩阵必须是 [[0,1],[-1,0]]。Ng 和 Np 分别指定焦点量和周期系数的最高指数。Ng=0 关闭焦点量；Ng>=2；Np=0 只保留周期常数。状态变量和待保留的符号参数请保持未赋值，精确计算建议使用整数或有理数。"),
        Code("X := (1+b*x)*(y+a*x*y):\nY := (1+b*x)*(-x+a*y^2):\nNg := 5:\nNp := 4:"),

        Heading("3. 计算并显示 G、H、g、p"),
        Code("ans := FocusPeriod:-Compute(X,Y,x,y,Ng,Np):\nFocusPeriod:-Show(ans);"),
        Text("本中心示例的预期结果：G[2]=a*sin(theta)，G[3]=a*b*cos(theta)*sin(theta)，H[0]=-1，H[1]=-b*cos(theta)，H[2]=0；g[2] 到 g[5] 全部为 0；p[2]=Pi*b*(b-a)，p[3]=-2*Pi*a*b*(b-a)。"),

        Heading("4. 查看中间量与周期展开"),
        Text("G、H、A、R 是按当前半径 r 展开的系数；u、du、B、p 是按初始半径 rho 展开的系数。A 是角速度倒数，B 是代入径向解后带负号的时间密度；p[n] 积分 B[n]，不能直接积分 A[n]。所有量都保存在 ans 中。theta、r、rho 是模块私有符号，代入前须从 ans 获取。u[n] 是表达式，初值通过 eval 查询。完整变量表、下标范围及批量查询见 README 第6节和 examples/03_inspect.mpl。"),
        Code("th := ans[\"theta\"]:\nrh := ans[\"rho\"]:\nans[\"A\"][0];\nans[\"A\"][1];\nans[\"R\"][2];\nans[\"u\"][2];\nans[\"du\"][2];\nsimplify(eval(ans[\"u\"][2],th=0));\nans[\"p\"][2];\nans[\"P_series\"];\neval(ans[\"P_series\"],{a=1,b=2});"),
        Text("最后一项应为 2*Pi+2*Pi*rh^2-4*Pi*rh^3+(15*Pi/2)*rh^4。u[2](0) 应为 0。"),

        Heading("5. 焦点示例：核对 g 的方向符号"),
        Text("此系统在 a 非零时为焦点；周期栏是形式时间展开。采用本工作表的角方向，g[3]=-a，g[5]=a*b+3*Pi*a^2。不要把有限个焦点量为零当成中心证明。"),
        Code("Xf := y+a*x*(x^2+y^2)+b*y*(x^2+y^2):\nYf := -x+a*y*(x^2+y^2)-b*x*(x^2+y^2):\nfocus := FocusPeriod:-Compute(Xf,Yf,x,y,5,4):\nFocusPeriod:-Show(focus);\nfocus[\"g\"][3];\nfocus[\"g\"][5];\nfocus[\"p\"][4];"),
        Text("可单独指定输出阶数，例如只算焦点量：FocusPeriod:-Compute(Xf,Yf,x,y,7,0)；只算中心的周期系数：FocusPeriod:-Compute(X,Y,x,y,0,6)。"),
        Text("参考算法：Zhang、Hou、Zeng (2000)，Computers & Mathematics with Applications 40，771–782，第 2 节与附录。本模块推广到任意多项式次数，并保留用户约定的 H[0]=-1 与 g 的方向。")
    );

    mwtext := :-Worksheet:-ToString(worksheet,format="mw");
    # 独立工作表必须定义 Maple Input 字体为 executable=true。
    # Layout 构造器适合插入已有GUI工作表；这里补齐独立文件所需样式。
    # 属性沿用 Maple 官方自带 examples/WorksheetPackage.mw 的标准样式。
    styles := cat(
        "<Version major=\"2025\" minor=\"0\"/>",
        "<View-Properties presentation=\"false\"/>",
        "<Styles>",
        "<Font name=\"Maple Input\" executable=\"true\" family=\"Courier New\" size=\"12\" foreground=\"[120,0,14]\" bold=\"true\"/>",
        "<Font name=\"Text\" executable=\"false\" family=\"Times New Roman\" size=\"12\" foreground=\"[0,0,0]\"/>",
        "<Font name=\"Heading 1\" executable=\"false\" family=\"Times New Roman\" size=\"20\" foreground=\"[0,0,0]\" bold=\"true\"/>",
        "<Font name=\"Heading 2\" executable=\"false\" family=\"Times New Roman\" size=\"16\" foreground=\"[0,0,0]\" bold=\"true\"/>",
        "<Layout name=\"Normal\" alignment=\"left\" linebreak=\"space\" spaceabove=\"0\" spacebelow=\"0\"/>",
        "<Layout name=\"Heading 1\" alignment=\"left\" linebreak=\"space\" spaceabove=\"8\" spacebelow=\"4\"/>",
        "<Layout name=\"Heading 2\" alignment=\"left\" linebreak=\"space\" spaceabove=\"8\" spacebelow=\"2\"/>",
        "</Styles>");
    mwtext := StringTools:-Substitute(mwtext,"<Worksheet>",cat("<Worksheet>",styles));
    FileTools:-Text:-WriteFile(outputFile,mwtext);

    # 使用 Maple 自己的解析器重读，及时发现生成文件格式错误。
    :-Worksheet:-ReadFile(outputFile,format="mw");
    printf("Native worksheet generated: %s\n",outputFile);
    return NULL;
end proc:

BuildWorksheet("examples/FocusPeriod_Worksheet.mw"):
