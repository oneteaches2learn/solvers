classdef GalerkinPennes2d_solver < GalerkinParabolic2d_solver

	properties
		%domain
		%time
		%coefficients
		%uInit
		%f
		%solution
	end

	methods
		function self = GalerkinPennes2d_solver(dom,time,cofs,uInit,f)
			
			% call superclass constructor
			self@GalerkinParabolic2d_solver(dom,time,cofs,uInit,f);

		end

		function self = initializeTensors(self,t)
			
			% call superclass method
			self = initializeTensors@GalerkinParabolic2d_solver(self);

			% create fields for tensor storage
			self.tensors.M_r = [];

			% check which tensors are time-varying
			self.tensors.timeVarying.M_r = Coefficients.isTimeVarying(self.coefficients.r);

		end

		function self = initializeVectors(self)

			% call superclass method
			self = initializeVectors@GalerkinParabolic2d_solver(self);

			% create fields for vector storage
			self.vectors.r_times_uStar = [];

			% check which vectors are time-varying
			self.vectors.timeVarying.r_times_uStar = ...
				(Coefficients.isTimeVarying(self.coefficients.uStar) || ...
				Coefficients.isTimeVarying(self.coefficients.r));
			self.vectors.timeVarying.r_times_uStar

		end

		function self = assembleTensors(self)

			% call superclass method
			self = assembleTensors@GalerkinParabolic2d_solver(self);

			% if first timestep, create tensors
			if self.isFirstTimeStep == 1
				self.tensors.M_r = self.assembleMassMatrix(self.coefficients.r);

			% else, update tensors as needed
			else
				% update M_r
				if self.tensors.timeVarying.M_r == 1
					fprintf('   Updating Tensor\n')
					cof = self.coefficients.r;
					self.tensors.M_r = self.assembleMassMatrix(cof);
				end
			end

		end

		function self = assembleVectors(self)

			% call superclass method
			self = assembleVectors@GalerkinParabolic2d_solver(self);

			% if first timestep, assemble r_times_uStar vector
			if self.isFirstTimeStep == 1
				self.vectors.r_times_uStar = self.compute_r_times_uStar;

			% else, update r_times_uStar vector as needed
			else
				if self.vectors.timeVarying.r_times_uStar == 1
					fprintf('    computing uStar\n')
					self.vectors.r_times_uStar = self.compute_r_times_uStar;
				end

			end
		end

		function b = compute_r_times_uStar(self)

			% store variables
			nNodes    = size(self.domain.Mesh.Nodes,2);
			nElem3    = size(self.domain.Mesh.Elements,2);
			coords    = self.domain.Mesh.Nodes';
			elements3 = self.domain.Mesh.Elements';
			uStar     = self.coefficients.uStar;
			r         = self.coefficients.r;

			% check fxn inputs
			if Coefficients.isTimeVarying(uStar) == 0
				uStar = @(x1,x2,t)(uStar(x1,x2));
			end

			if Coefficients.isTimeVarying(r) == 0
				r = @(x1,x2,t)(r(x1,x2));
			end

			% initialize storage
			b = sparse(nNodes,1);

			% compute volume forces
			for j = 1:nElem3
				elementInd    = elements3(j,:);
				elementCoord  = coords(elementInd,:);
				b(elementInd) = b(elementInd) + ...
					det([1,1,1; elementCoord']) * ...
					r(sum(elementCoord(:,1))/3,sum(elementCoord(:,2))/3,self.t) * ...
					uStar(sum(elementCoord(:,1))/3,sum(elementCoord(:,2))/3,self.t) / 6;

			end
		end

		function [S,b] = finalAssembly(self,U_prev)

			% store variables
			tensors = self.tensors;
			vectors = self.vectors;

			% assemble LHS
			S = self.time.dt * (tensors.A + tensors.M_r + tensors.E) + tensors.M_p;

			% assemble RHS
			b = self.time.dt * (vectors.b_vol - vectors.b_neu + vectors.b_rob + ...
				 		 vectors.r_times_uStar) - S * vectors.U_D + tensors.M_p_prev * U_prev;

		end

	end
end
