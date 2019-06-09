#!/usr/bin/env luajit
-- todo paths for non-windows
for _,info in ipairs{
	--{name='gauss_seidel', args='selfGravPoissonSolver=poisson'},
	{name='conjgrad', args='selfGravPoissonSolver=poisson_krylov selfGravLinearSolver=conjgrad'},
	{name='conjres', args='selfGravPoissonSolver=poisson_krylov selfGravLinearSolver=conjres'},
	{name='gmres', args='selfGravPoissonSolver=poisson_krylov selfGravLinearSolver=gmres'},
} do
	os.execute([[cd ..\.. && luajit run.lua sys=console selfGravVerbose exitTime=0 selfGravInitPotentialPositive ]]..info.args..[[ 2> tests\poisson-convergence\init-]]..info.name..[[-pos.txt"]])
	os.execute([[cd ..\.. && luajit run.lua sys=console selfGravVerbose exitTime=0 ]]..info.args..[[ 2> tests\poisson-convergence\init-]]..info.name..[[-neg.txt"]])
end
