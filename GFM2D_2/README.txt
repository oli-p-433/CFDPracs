Initial conditions, BCs, solver, slope limiter and riemann solver can be changed within fvm.C (and CFL)

solver options are MUSCL, godunov, SLIC, FORCE, Richt, 

Comment out all irrelevant initial conditions to select test case (very janky, sorry)

Save fvm.C and run ./allrun ; results will be stored in the folder t

to view results, use command:

./plot {variable = rho1, rho2, vx1, vx2, vy1, vy2, p1, p2} {folder} {timestep}

eg ./plot rho1 t 0.100
(if in doubt check folder for available timesteps)
