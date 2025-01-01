classdef NewtonGalerkinSolver2d_elliptic < GalerkinSolver2d_elliptic & NewtonGalerkinSolver2d

	properties
		U
	end

	methods
		function self = NewtonGalerkinSolver2d_elliptic(dom,auxfun)
			
			% call superclass constructor
			self@GalerkinSolver2d_elliptic(dom,auxfun);
			self@NewtonGalerkinSolver2d();

		end

		function self = solve(self)

			% Initial guess
			self.U = ones(self.domain.mesh.nNodes,1);

			% Initialisation
			dirichlet = self.domain.boundary.D_nodes;
			FreeNodes = self.domain.boundary.freeNodes;
			self = self.assembleTensors;
			self = self.assembleVectors;

			% Update initial guess with Dirichlet BC values
			self = self.assembleBCs;
			self.U(unique(dirichlet)) = self.vectors.U_D(unique(dirichlet));

			% temporary
			%self.U(unique(dirichlet)) = 0;

			% Newton-Galerkin iteration
			for i = 1:100
			
				% Assembly
				self = self.assembleTensors;
				self = self.assembleVectors;
				self = self.assembleBCs;
				[DJ,J] = self.finalAssembly(self.U);
							
				% Dirichlet conditions
				W = zeros(self.domain.mesh.nNodes,1);
				W(unique(dirichlet)) = 0;
				
				% Solving one Newton step
				W(FreeNodes) = DJ(FreeNodes,FreeNodes) \ J(FreeNodes);
				self.U = self.U - W;

				% check convegence
				if norm(W) < 10^(-10)
					fprintf(' %d iterations,',i)
					break
				end

			end

			% store result
			self.solution = self.U;

		end


		function [DJ,J] = finalAssembly(self,U)

			% note: placeholder function. Actually handled by specific subclasses.
			% ...

		end

		function self = assembleTensors(self) 

			% call superclass method
			self = self.assembleTensors@GalerkinSolver2d_elliptic;

			% assemble additional tensors
			self.tensors.M_dr = self.assembleMassMatrix(self.coefficients.dr_du);
			self.tensors.M_dneu = self.computeNonlinearNeumannContribution;
			self.tensors.M_drob = self.computeNonlinearRobinContribution;

		end

		function self = assembleVectors(self) 

			% call superclass method
			self = self.assembleVectors@GalerkinSolver2d_elliptic;

			% assemble additional tensors
			self.vectors.M_r = self.computeVolumeForces(self.coefficients.r);

		end

		function self = assembleBCs(self)

			% note: as of Oct 11, 2024 periodic BCs don't work, so they are muted
			%self = self.computePeriodicBCs;
			self.vectors.U_D   = self.computeDirichletBCs;
			self.vectors.b_neu = self.computeNeumannBCs;
			[self.tensors.M_rob,self.vectors.b_rob] = self.computeRobinBCs;
			%[self.tensors.M_rob,self.tensors.M_dyn,self.vectors.b_dyn] = self.computeDynamicBCs;

		end

	end
end



		% NOTE: This is all old code that has been muted out at some point. I
		% don't know if it is worth keeping any of this. I think it's all just
		% deprecated and useless, or copies of functions that appear elsewhere.
		% So consider the following a graveyard.  
		%{
		function A = assembleStiffnessMatrix(self,U)

			% Initialisation
			nodes = self.domain.mesh.nodes;
			elems = self.domain.mesh.elements;

			A = sparse(self.domain.mesh.nNodes,self.domain.mesh.nNodes);
			for j = 1:self.domain.mesh.nElems
				A(elems(j,:),elems(j,:)) = A(elems(j,:),elems(j,:)) ...
					+ self.localdj(nodes(elems(j,:),:),U(elems(j,:)));
			end

		end
		%}

		%{
		function self = assembleVectors(self,U)

			% Initialisation
			nodes = self.domain.mesh.nodes;
			elems = self.domain.mesh.elements;
			[temp,b_cubic] = self.domain.nodalQuadrature(U.^3);

			b = sparse(size(nodes,1),1);
			for j = 1:self.domain.mesh.nElems
				b(elems(j,:)) = b(elems(j,:)) ...
					+ self.localj(nodes(elems(j,:),:),U(elems(j,:)),b_cubic(j));
			end

			self.vectors.b_diffeq = b;
			self.vectors.b_load = self.computeVolumeForces;

		end

		function self = assembleBCs(self)

			%self = self.computePeriodicBCs;
			self.vectors.U_D   = self.computeDirichletBCs;
			self.vectors.b_neu = self.computeNeumannBCs;
			[self.tensors.E,self.vectors.b_rob] = self.computeRobinBCs;

		end

		function M = localdj(self,vertices,U_loc)

			% compute local contributions to DJ
			A_loc = self.localStiffnessMatrix(vertices);
			M_loc = self.localMassMatrix(vertices);
			b_Dcubic = self.localDCubic_direct(U_loc);
			b_Dcubic = self.localDCubic_quad(U_loc);

			% assemble local contributions to DJ
			M = A_loc - M_loc + b_Dcubic;
			Area = det([ones(1,3);vertices']) / 2;
			M = M * Area;
			
		end
		
		function j_loc = localj(self,vertices,U_loc,U_cubed)

			% compute local contributions to J
			A_loc = self.localStiffnessMatrix(vertices);
			b_nonlinear = self.localQuadrature(U_loc.^3 - U_loc);

			% assemble local contributions to J
			Area = det([ones(1,3); vertices']) / 2;
			j_loc = A_loc * U_loc + b_nonlinear;
			j_loc = Area * j_loc;

		end

		function A_loc = localStiffnessMatrix(self,vertices)

			% compute local stiffness matrix
			G = [ones(1,3);vertices'] \ [zeros(1,2);eye(2)];
			A_loc = G * G';

		end

		function M_loc = localMassMatrix(self,vertices)

			% compute local mass matrix
			M_loc = [2,1,1; 1,2,1; 1,1,2] / 12;

		end
		
		function b_cubic = localCubic_quad(self,U)

			M_loc = [2,1,1; 1,2,1; 1,1,2] / 12;
			b_cubic = M_loc * [U(1)^3; U(2)^3; U(3)^3];

		end

		function b = localQuadrature(self,U)

			M_loc = [2,1,1; 1,2,1; 1,1,2] / 12;
			b = M_loc * U;
			
		end

		function b_cubic = localCubic_direct(self,U)

			% direct computation of cubic integral
			b_cubic = ...
			[4 * U(1)^3 + U(2)^3 + U(3)^3 + 3 * U(1)^2 * (U(2) + U(3)) + 2 * U(1) ...
				* (U(2)^2 + U(3)^2) + U(2) * U(3) * (U(2) + U(3)) + 2 * U(1) * U(2) * U(3);
			 4 * U(2)^3 + U(1)^3 + U(3)^3 + 3 * U(2)^2 * (U(1) + U(3)) + 2 * U(2) ...
				* (U(1)^2 + U(3)^2) + U(1) * U(3) * (U(1) + U(3)) + 2 * U(1) * U(2) * U(3);
			 4 * U(3)^3 + U(2)^3 + U(1)^3 + 3 * U(3)^2 * (U(2) + U(1)) + 2 * U(3) ...
				* (U(2)^2 + U(1)^2) + U(2) * U(1) * (U(2) + U(1)) + 2 * U(1) * U(2) * U(3)] / 60;

		end


		function b_Dcubic = localDCubic_direct(self,U)

			% direct computation of cubic derivative integral
			b_Dcubic = ...
			[12 * U(1)^2 + 2 * (U(2)^2 + U(3)^2 + U(2) * U(3)) + 6 * U(1) * (U(2) + U(3)), ...
				 3 * (U(1)^2 + U(2)^2) + U(3)^2 + 4 * U(1) * U(2) + 2 * U(3) * (U(1) + U(2)), ...
			  	 3 * (U(1)^2 + U(3)^2) + U(2)^2 + 4 * U(1) * U(3) + 2 * U(2) * (U(1) + U(3));
			  3 * (U(1)^2 + U(2)^2) + U(3)^2 + 4 * U(1) * U(2) + 2 * U(3) * (U(1) + U(2)), ...
				12 * U(2)^2 + 2 * (U(1)^2 + U(3)^2 + U(1)*U(3)) + 6 * U(2) * (U(1) + U(3)), ...
				 3 * (U(2)^2 + U(3)^2) + U(1)^2 + 4 * U(2) * U(3) + 2 * U(1) * (U(2) + U(3));
		 	  3 * (U(1)^2 + U(3)^2) + U(2)^2 + 4 * U(1) * U(3) + 2 * U(2) * (U(1) + U(3)), ...
				 3 * (U(2)^2 + U(3)^2) + U(1)^2 + 4 * U(2) * U(3) + 2 * U(1) * (U(2) + U(3)), ...
				12 * U(3)^2 + 2 * (U(1)^2 + U(2)^2 + U(1) * U(2)) + 6 * U(3) * (U(1) + U(2))] / 60;

		end

		function b_Dcubic = localDCubic_quad(self,U)

			M_loc = [2,1,1; 1,2,1; 1,1,2] / 12;
			b_Dcubic = 3 * (U) * (U)' .* M_loc;

		end

		function b = computeVolumeForces(self,U)

			% store variables
			f = self.f;

			% check coefficient variables
			if Coefficients.isTimeVarying(f) == 0 && Coefficients.isNonlinear(f) == 0
				f = @(x1,x2,t,u)(f(x1,x2));
				u = 0;
                t = 0;
			elseif Coefficients.isTimeVarying(f) == 1 && Coefficients.isNonlinear(f) == 0
				f = @(x1,x2,t,u)(f(x1,x2,t));
				u = 0;
			elseif Coefficients.isTimeVarying(f) == 0 && Coefficients.isNonlinear(f) == 1
				f = @(x1,x2,t,u)(f(x1,x2,u));
                t = 0;
            else
                t = self.t;
			end

			% unpack variables
			nNodes    = self.domain.mesh.nNodes;
			nElem3    = self.domain.mesh.nElems;
			coords    = self.domain.mesh.nodes;
			elements3 = self.domain.mesh.elements;

			% initialize storage
			b = sparse(nNodes,1);

			% compute volume forces
			for j = self.domain.mesh.effectiveElems
				elementInd    = elements3(j,:);
				elementCoord  = coords(elementInd,:);
				b(elementInd) = b(elementInd) + ...
					self.domain.mesh.areas(j) * ...
					f(sum(elementCoord(:,1))/3,sum(elementCoord(:,2))/3,t,U(elementInd)) / 3;
			end
		end
		%}

