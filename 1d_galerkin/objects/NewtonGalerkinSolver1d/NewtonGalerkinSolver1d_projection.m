classdef NewtonGalerkinSolver1d_projection < NewtonGalerkinSolver1d_elliptic
% NEWTONGALERKINSOLVER1D_PROJECTION finds the solution W to the elliptic projection
%
%  (k nabla (u - W), nabla z)_Om + lambda (u - W,z)_Om - (a(u,v) - a(W,v),z)_pOm = 0
%
% for all test functions z, where (u,v) are known functions, where lambda > 0,
% and where a is a nonlinear boundary term. 
%
% This projection operator could be thought of as the weak form of the PDE
%
% -nabla . (k nabla W) + lambda W = -nabla . (k nabla u) + lambda u, in Omega
%                  -k nabla W . n = a(W,v),                          on dOmega
%
% Therefore, this PDE is suitably solved using the NewtonGalerkinSolver1d class.
% In particular, one could write the PDE as
% 
% -nabla . (k nabla W) + lambda W = f,        in Omega
%                  -k nabla W . n = a(W,v),   on dOmega
%
% with f = -nabla . (k nabla u) + lambda u. The NewtonGalerkinSolver1d class can
% already handle nonlinear boundary condtions of the form g = g(W,v). So the
% only additional feature that needs to be added is to handle the source term f.
%
% To handle f, we discretize the entire PDE using the Galerkin method, leading
% to the system in matrix form
% 
% 	(A + lambda M) W + M_neu(W,v) = (A + lambda M) U + b_neu(U,v)           (*)
%
% where A is the stiffness matrix, M is the mass matrix, M_neu is the nonlinear
% boundary term, and b_neu is the load vector corresponding to the boundary
% term. U is a nodal interpolation of u. The resulting right-hand side replaces
% the usual b_vol vector. And then the Newton--Galerkin method can proceed as
% usual. Notice that the NewtonGalerkinSolver1d class can already compute all of
% the necessary matrices and vectors. So we just need to modify the
% finalAssembly method.
%
%  Author: Tyler Fara                             Date: November 16, 2025
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% NOTES:
%   (1) The usual auxfun input should have auxfun.u_fixed and auxfun.v_fixed in
%   place of the usual auxfun.f. As usual, these should be function handles of
%   space and time. 
%
%   (2) Also, auxfun should have auxfun.cofs.lambda with lambda > 0.
%
%   (3) The known functions u and v are saved as u_fixed and v_fixed. The
%   subscript "fixed" indicates that these functions are held constant during
%   the projection. (That is, in this class, the term "fixed" is not a reference
%   to anything having been "broken".)
%
%   (4) The superclass NewtonGalerkinSolver1d uses property U to store the
%   current solution guess. Therefore, when you see self.U in this class, it
%   actually refers to the current guess W.
%
%   (5) Robin boundary conditions have been omitted for now. In other classes, I
%   have kept commented out versions of Robin BCs for reference. This is because
%   Robin conditions have been successfully implemented in the GalerkinSolver2d
%   family. Since the GalerkinSolver1d classes were made by modifying their 2d
%   counterparts, I am reasonably sure that the implementation of Robin
%   conditions is the other GalerkinSolver1d classes is correct; it's just not
%   complete. However, in this class, I'm not actually sure how Robin conditions
%   would be implemented. They might require a totally different implementation,
%   or maybe they are impossible. At any rate, a Robin condition could be
%   instituted by using the nonlinear Neumann condition anyway. So, to avoid
%   possible errors, I have simply removed all references to Robin conditions.
%   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



	properties
	end

	methods
		function self = NewtonGalerkinSolver1d_projection(dom,auxfun)
			
			% call superclass constructor
			self@NewtonGalerkinSolver1d_elliptic(dom,auxfun);

			% evaluate known functions on mesh nodes
			nodes = self.domain.mesh.nodes;
			U_fixed = self.coefficients.u_fixed(nodes);
			V_fixed = self.coefficients.v_fixed(0);

			% set initial guess to (U_fixed,V_fixed)
			self.U = U_fixed;
			self.V = V_fixed;

			% store RHS
			self.vectors.b_neuFixed = self.computeNeumannBCs(U_fixed,V_fixed);
			self = self.assembleTensors;
			self = self.assembleVectors;
			self.vectors.b_vol = self.computeRHS(U_fixed,V_fixed);

			% temporary: overwrite b_vol with the exact value, computed by hand
			%self.vectors.b_vol_estimate = self.vectors.b_vol;
			%self.vectors.b_vol = [-1/2; -1; 3/2] + 1 * [1/96; 7/48; 17/96] - [0; 0; 1];

			% solve
			%self = self.solve;

		end

		function [DJ,J] = finalAssembly(self,U_tilde)
		% FINALASSEMBLY performs the final assembly of the system matrices/vectors.
		%
		% NOTES:
		%   (1) This is actually the exact same code as the
		%   NewtonGalerkinSolver1d_poisson subclass. The difference is that the
		%   b_vol term here represents the RHS of (*) from the class header,
		%   instead of the usual load vector. But with that change, the assembly
		%   is otherwise exactly the same. (That is, minus the fact that we make
		%   no attempt to implement Robin BCs here.)
		

		% NOTE: I am changing the sign on the Neumann conditions to match the
		% sign conentions of my paper. If you see "+" signs by the b_neu and
		% M_dneu terms, then this is after the change. It originally had minus
		% signs.  

			% store variables
            tensors = self.tensors;
            vectors = self.vectors;

			% assemble Linear Tensor
			S = tensors.A;

			% assemble nonlinear contributions to J
			b_nonlinear = self.vectors.M_r;

			% assemble Load Vector
			b = (vectors.b_vol + vectors.b_neu) - S * vectors.U_D;

			% assemble J
            J = S * U_tilde + b_nonlinear - b;
                
			% assemble DJ
			%DJ = tensors.A + self.coefficients.lambda * tensors.M_dr + tensors.M_dneu;
			DJ = tensors.A + tensors.M_dr - tensors.M_dneu;

		end


		function self = assembleVectors(self)
		% ASSEMBLEVECTORS assembles the RHS vectors and overwrites the superclass method. 
		% 
		% In particular, b_vol now represents the right hand side of (*) in the
		% description of this class, instead of the usual load vector. The
		% vector M_r is still assembled as usual.
		%
		% NOTE: b_vol is data. You actually do not want to recompute it every
		% time. Indeed, it is now set by the constructor, which is good because
		% it means you can overwrite it before calling solve.
		
			%self.vectors.b_vol = self.computeRHS(self.u_fixed,self.v_fixed);
			self.vectors.M_r = self.computeVolumeForces(self.coefficients.r);

		end

		function self = assembleBCs(self)
		% ASSEMBLEBCS assembles the boundary condition contributions and overwrites
		% the superclass method. 
		%
		% NOTES:
		%   (1) The method computeNeumannBCs is called twice: once to compute
		%   the boundary vector b_neu corresponding to the current guess W
		%   (which is as usual), and once to compute the vector b_neuFixed.
		%
		%   (2) In the NewtonGalerkinSolver1d superclass, property U is the
		%   current guess. Therefore, when we pass self.U below, we are actually
		%   passing the current guess for solution W.
		
			% store variables
			V_fixed = self.coefficients.v_fixed(0);

			self.vectors.U_D   = self.computeDirichletBCs;
			self.vectors.b_neu = self.computeNeumannBCs(self.U,V_fixed);
			%[self.tensors.M_rob,self.vectors.b_rob] = self.computeRobinBCs;

		end

		function b_vol = computeRHS(self,U_fixed,V_fixed)
		% COMPUTERHS computes the right-hand side vector corresponding to the
		% projection problem description in the class header.

			% load variables
			tensors = self.tensors;
			vectors = self.vectors;

			% input handling
			if nargin < 2
				U_fixed = self.coefficients.u_fixed(self.domain.mesh.nodes);
				V_fixed = self.coefficients.v_fixed(0);
			end

			% compute right-hand side vector
			%b_vol = (tensors.A + lambda * tensors.M_r) * U_fixed - vectors.b_neuFixed;
			b_vol = (tensors.A + tensors.M_r) * U_fixed - vectors.b_neuFixed;

		end

		function b_neu = computeNeumannBCs(self,U,V)
		% COMPUTENEUMANNBCS computes the Neumann boundary conditions contribution to the
		% right-hand side vector. This method overwrites the superclass method
		% and allows to pass U and V as arguments.
		%
		% NOTES:
		%   (1) In the superclass (i.e. GalerkinSolver1d), it is assumed that U
		%   and V would correspond to the current guess for the solution values.
		%   However, in this projection problem, we actually need to compute two
		%   different Neumann boundary vectors: 
		%   
		% 	 (a) one for the current guess W, used by the Newton iteration, and
		% 	 (b) one for the known function u, used to compute the right-hand side.
		%
		%   By allowing U and V to be passed as arguments, we can compute both
		%   of these with the same code.
		
            % unpack variables
            dom    = self.domain;
            nNodes = dom.mesh.nNodes;
            coords = dom.mesh.nodes(:);   % ensure column

            % initialize storage
            b_neu = sparse(nNodes, 1);

            % loop over boundary edges
            for i = 1:dom.boundary.nEdges

                edge = dom.boundary.edges(i);

                % only process Neumann edges
                if edge.boundaryType == 'N'

                    % nodes on this Neumann boundary edge
                    bNodes_i = edge.nodes;          % indices into global node list
                    bCond    = edge.boundaryCondition;

                    % normalize boundary condition into f(x,t,u,v) form
                    [bCond, t] = self.checkVariables(bCond);

                    % coordinates and solution values at these boundary nodes
                    x_b   = coords(bNodes_i);
                    u_b   = U(bNodes_i);

                    % evaluate Neumann flux at boundary nodes
                    q = bCond(x_b, t, u_b, V);      % could be scalar or vector

                    if numel(q) == 1
                        q = q * ones(numel(bNodes_i), 1);  % broadcast scalar
                    end

                    % accumulate into global RHS
                    b_neu(bNodes_i) = b_neu(bNodes_i) + q;
                end
            end
        end



	end
end
