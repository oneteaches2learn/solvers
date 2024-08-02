(1) To demo, in MATLAB, run: 

  packageMMS_demo

This runs an MMS (method of manufactured solutions) test of the 2d Galerkin
Pennes equation solver. The solver will create an mms object containing four
solved trials with successively refined grid, demonstrating order 2
convergence. 

(2) After running package_demo, try:

  mms.plot
  mms.problems{4}.domain.plot
  mms.problems{4}.animate

These will show, respectively: a convergence plot for the MMS test, the mesh
for the most refined problem, and an animation of the solution for the most
refined problem.

(3) Open solvers > packageMMS_demo. Try changing parameters under "USER
INPUTS". For example, try changing the desired solution to be manufactured,
uTrue. Input any desired function using the variables x(1), x(2), t. 

(4) To demo the solver itself, run:

  packageSolve_demo

Then try:

  prob.animate
  prob.domain.plot
  prob.getSolutionIntegral

These will show, respectively: an animation of the solution, the mesh, and a
vector of the L2 norm of the solution at each timestep. 

(5) Open solvers > packageSolve_demo. Try changing parameters. Different
boundary conditions can be specified on different edges of the domain using 'D'
for dirichlet, 'N' for neumann, and 'R' for Robin.

(6) Inspect the base class GalerkinParabolic2d_solver, found under: 

  solvers > 2d_galerkin > 2d_parabolic_linear > objects > GalerkinParabolic2d_solver.m

Specific problems (e.g. the heat equation, the reaction diffusion equation, the
pennes equation) are solved by subclasses of GalerkinParabolic2d_solver. For
example, the pennes equation is solved by GalerkinPennes2d_solver, found under:

  solvers > 2d_galerkin > 2d_pennes_linear > objects > GalerkinPennes2d_solver.m
