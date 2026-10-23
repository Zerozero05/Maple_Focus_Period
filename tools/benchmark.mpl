# 在仓库根目录用cmaple -q tools/benchmark.mpl运行；仅测量，不写文件。
# 每种设置先预热，再取三次CPU/墙钟时间中位数。结果依赖本机与具体模型。
restart:
read "src/focus_period.mpl":
Benchmark := proc(label::string,X,Y,Ng::nonnegint,Np::nonnegint)
    local reference, sample, variant, repetition, samples_cpu, samples_wall,
          started_cpu, started_wall, n;
    reference := FocusPeriod:-Compute(X,Y,x,y,Ng,Np,method="quotient",
                                     optimize=false,factor_orbits=false);
    for variant in ["quotient_dense","direct_dense","direct_default"] do
        samples_cpu := []; samples_wall := [];
        for repetition from 0 to 3 do
            started_cpu := time(); started_wall := time[real]();
            if variant="quotient_dense" then
                sample := FocusPeriod:-Compute(X,Y,x,y,Ng,Np,method="quotient",
                                              optimize=false,factor_orbits=false);
            elif variant="direct_dense" then
                sample := FocusPeriod:-Compute(X,Y,x,y,Ng,Np,
                                              optimize=false,factor_orbits=false);
            else
                sample := FocusPeriod:-Compute(X,Y,x,y,Ng,Np);
            end if;
            if repetition>0 then
                samples_cpu := [op(samples_cpu),time()-started_cpu];
                samples_wall := [op(samples_wall),time[real]()-started_wall];
            end if;
            for n from 2 to Ng do
                if not evalb(simplify(sample["g"][n]-reference["g"][n])=0) then
                    error "Benchmark g mismatch: %1, %2, order %3",label,variant,n;
                end if;
            end do;
            for n from 0 to Np do
                if not evalb(simplify(sample["p"][n]-reference["p"][n])=0) then
                    error "Benchmark p mismatch: %1, %2, order %3",label,variant,n;
                end if;
            end do;
        end do;
        printf("%s,%s,%d,%d,3,%.6f,%.6f,verified\n",label,variant,Ng,Np,
               sort(samples_cpu)[2],sort(samples_wall)[2]);
    end do;
end proc:
printf("model,configuration,Ng,Np,repeats,median_cpu_s,median_wall_s,g_p\n"):
Benchmark("radial_cubic",y+a*x*(x^2+y^2)+b*y*(x^2+y^2),
          -x+a*y*(x^2+y^2)-b*x*(x^2+y^2),11,10):
Benchmark("factored_center",(1+b*x)*(y+a*x*y),(1+b*x)*(-x+a*y^2),9,8):
Benchmark("mixed_quadratic_cubic",y+x^2+x*y+2*y^2+x^3,-x+x*y+y^2+x*y^2,7,6):
printf("BENCHMARK PASS: all measured g/p agree with the quotient reference.\n"):
