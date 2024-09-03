classdef (Abstract) GalerkinSolver2d_elliptic < GalerkinSolver2d

	properties
	end

	methods
		function self = GalerkinSolver2d_elliptic(dom,cofs)
			
			% call superclass constructor
			self@GalerkinSolver2d(dom,cofs);

			% calculate solution
			% ...

		end

		function self = solve(self)

			% load variables
			FreeNodes = self.domain.freeNodes;

			% assemble problem
			self = self.assembleTensors;
			self = self.assembleVectors;
			self = self.assembleBCs;
			[S,b] = self.finalAssembly;

			% solve and store solution
			v = sparse(self.domain.mesh.nNodes,1);
			v(FreeNodes) = S(FreeNodes,FreeNodes) \ b(FreeNodes);
			self.solution = v + self.vectors.U_D;

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

		function [S,b] = finalAssembly(self)

			% NOTE: placeholder function. Actually handled by specific subclasses.
			...

		end

	end
end
