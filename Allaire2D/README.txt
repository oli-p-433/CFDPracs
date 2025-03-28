Initial conditions, BCs and solver can be changed within fvm.C (and CFL)

Solver options are HLLCGodunov, MUSCL, SLIC.
set bool prim = true; for primitive MUSCL, false for conservative MUSCL

Comment out all irrelevant initial conditions to select test case (very janky, sorry)

Save fvm.C and run ./allrun ; results will be stored in the folder t

to view results, use command:

./plot {variable} {folder} {timestep}

for all variables apart from rho; for rho, use

./plotRho {folder} {timestep}
