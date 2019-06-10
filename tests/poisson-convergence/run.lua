#!/usr/bin/env luajit
-- todo paths for non-windows

local function exec(cmd)
	print(cmd)
	return os.execute(cmd)
end

local commonArgs = 'sys=console eqn=euler solver=roe "initState=self-gravitation test 1" dim=2 selfGravVerbose exitTime=0 selfGravPoissonMaxIter=10000'

for _,info in ipairs{
	{name='jacobi', args='selfGravPoissonSolver=jacobi'},
	{name='conjgrad', args='selfGravPoissonSolver=krylov selfGravLinearSolver=conjgrad'},
	{name='conjres', args='selfGravPoissonSolver=krylov selfGravLinearSolver=conjres'},
	{name='gmres', args='selfGravPoissonSolver=krylov selfGravLinearSolver=gmres'},
} do
	exec([[cd ..\.. && luajit run.lua ]]..commonArgs..' '..info.args..[[ selfGravInitPotential=+ 2> tests\poisson-convergence\init-]]..info.name..[[-pos.txt"]])
	exec([[cd ..\.. && luajit run.lua ]]..commonArgs..' '..info.args..[[ 2> tests\poisson-convergence\init-]]..info.name..[[-neg.txt"]])
end
