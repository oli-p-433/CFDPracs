Initial conditions, BCs and solver can be changed within fvm1D.C

Comment out all irrelevant initial conditions to select test case (very janky, sorry)

Save fvm1D.C and run ./allrun ; results will be stored in the folder t

to view results, use command:

./plot {variable} {folder} {timestep}

for all variables apart from rho; for rho, use

./plotRho {folder} {timestep}
