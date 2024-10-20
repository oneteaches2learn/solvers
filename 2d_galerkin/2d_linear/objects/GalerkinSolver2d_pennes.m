classdef GalerkinSolver2d_pennes < GalerkinSolver2d_parabolic

	methods
		function self = GalerkinSolver2d_pennes(dom,auxfun)
			
			% call superclass constructor
			self@GalerkinSolver2d_parabolic(dom,auxfun);

		end

		function self = initializeTensors(self,t)
			
			% call superclass method
			self = initializeTensors@GalerkinSolver2d_parabolic(self);

			% create fields for tensor storage
			self.tensors.M_r = [];

			% check which tensors are time-varying
			self.tensors.timeVarying.M_r = Coefficients.isTimeVarying(self.coefficients.r);

		end

		function self = initializeVectors(self)

			% call superclass method
			self = initializeVectors@GalerkinSolver2d_parabolic(self);

			% create fields for vector storage
			self.vectors.r_times_uStar = [];

			% check which vectors are time-varying
			self.vectors.timeVarying.r_times_uStar = ...
				(Coefficients.isTimeVarying(self.coefficients.uStar) || ...
				Coefficients.isTimeVarying(self.coefficients.r));

		end

		function self = assembleTensors(self)

			% call superclass method
			self = assembleTensors@GalerkinSolver2d_parabolic(self);

			% if first timestep, create tensors
			if self.isFirstTimeStep == 1
				self.tensors.M_r = self.assembleMassMatrix(self.coefficients.r);

			% else, update tensors as needed
			else
				% update M_r
				if self.tensors.timeVarying.M_r == 1
					cof = self.coefficients.r;
					self.tensors.M_r = self.assembleMassMatrix(cof);
				end
			end

		end

		function self = assembleVectors(self)

			% call superclass method
			self = assembleVectors@GalerkinSolver2d_parabolic(self);

			% if first timestep, assemble r_times_uStar vector
			if self.isFirstTimeStep == 1
				self.vectors.r_times_uStar = self.compute_r_times_uStar;

			% else, update r_times_uStar vector as needed
			else
				if self.vectors.timeVarying.r_times_uStar == 1
					self.vectors.r_times_uStar = self.compute_r_times_uStar;
				end

			end
		end

		function b = compute_r_times_uStar(self)

			% store variables
			uStar     = self.coefficients.uStar;
			r         = self.coefficients.r;

			% check fxn inputs
			if Coefficients.isTimeVarying(uStar) == 0
				uStar = @(x1,x2,t)(uStar(x1,x2));
			end

			if Coefficients.isTimeVarying(r) == 0
				r = @(x1,x2,t)(r(x1,x2));
			end

			% compute volume forces
			b = self.computeVolumeForces(uStar,r);

			%{
			% OLD CODE: loop-based version
			% store more variables
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
					det([1,1,1; elementCoord']) * ...
					r(sum(elementCoord(:,1))/3,sum(elementCoord(:,2))/3,self.t) * ...
					uStar(sum(elementCoord(:,1))/3,sum(elementCoord(:,2))/3,self.t) / 6;

			end
			%}

		end

		function [S,b] = finalAssembly(self)

			% store variables
			tensors = self.tensors;
			vectors = self.vectors;

			% assemble LHS
			S = self.domain.time.dt * (tensors.A + tensors.M_r + tensors.M_rob + tensors.M_dyn_u) + ...
							tensors.M_p + tensors.M_dyn_du;

			% assemble RHS
			b = self.domain.time.dt * (vectors.b_vol - vectors.b_neu + vectors.b_rob + ...
			 		 vectors.r_times_uStar + vectors.b_dyn) ...
					 - S * vectors.U_D + (tensors.M_p_prevTime + tensors.M_dyn_prevTime) * vectors.U_prevTime;

		end

	end
end
