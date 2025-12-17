classdef (Abstract) GalerkinSolver1d_elliptic < GalerkinSolver1d

	properties (Hidden)
	end

	methods
		function self = GalerkinSolver1d_elliptic(dom,auxfun)
			
			% call superclass constructor
			self@GalerkinSolver1d(dom,auxfun);

			if nargin == 2

				% calculate solution
				self.domain.time = [];
            	%self = self.solve;

			end

		end

		function self = solve(self)

			% assemble problem
			self = self.assembleTensors;	
			self = self.assembleVectors;
			self = self.assembleBCs;
			[S,b] = self.finalAssembly;
			%[S,b] = self.periodicCorrection(S,b);

			% load variables
			FreeNodes = self.domain.boundary.freeNodes;

			% solve and store solution
			v = sparse(self.domain.mesh.nNodes,1);
			v(FreeNodes) = S(FreeNodes,FreeNodes) \ b(FreeNodes);
			self.solution = v + self.vectors.U_D;

            %{
			% copy solution to periodic replica nodes
			P = self.domain.boundary.P_nodes;
			self.solution(P.replica.edge,:) = self.solution(P.free.edge,:);
			self.solution(P.replica.corner,:) = ...
							repmat(self.solution(P.free.corner,:),3,1);					
            %}

			% cleanup
			self = self.cleanup;

		end

		function self = assembleTensors(self) 

			self.tensors.A = self.assembleStiffnessMatrix;
			self.tensors.M_r = self.assembleMassMatrix(self.coefficients.r);

		end

		function self = assembleVectors(self)

			self.vectors.b_vol = self.computeVolumeForces();

		end

		function self = assembleBCs(self)

			%self = self.computePeriodicBCs;
			self.vectors.U_D   = self.computeDirichletBCs;
			self.vectors.b_neu = self.computeNeumannBCs;
			%[self.tensors.M_rob,self.vectors.b_rob] = self.computeRobinBCs;

		end

		function [S,b] = finalAssembly(self)

			% NOTE: placeholder function. Actually handled by specific subclasses.
			...

		end

	end
end
