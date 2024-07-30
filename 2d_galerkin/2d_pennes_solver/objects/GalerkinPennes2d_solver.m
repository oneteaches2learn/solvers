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
			
			% Call superclass method
			self = initializeTensors@GalerkinParabolic2d_solver(self);

			% create fields for tensor storage
			self.tensors.M_r = [];

			% check which tensors are time-varying
			self.tensors.timeVarying.M_r = Coefficients.checkTimeVarying(self.coefficients.r);

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
					cof = self.coefficients.r;
					self.tensors.M_r = self.assembleMassMatrix(cof);
				end
			end

		end

		function self = assembleVectors(self)

			% call superclass method
			self = assembleVectors@GalerkinParabolic2d_solver(self);

			% assemble uStar vector
			self.vectors.uStar = self.compute_uStar;

		end

		function [S,b] = finalAssembly(self,U_prev)

			% store variables
			tensors = self.tensors;
			vectors = self.vectors;

			% assemble LHS
			S = self.time.dt * (tensors.A + tensors.M_r + tensors.E) + tensors.M_p;

			% assemble RHS
			b = self.time.dt * (vectors.b_vol - vectors.b_neu + vectors.b_rob + ...
					vectors.uStar) - S * vectors.U_D + tensors.M_p_prev * U_prev;

		end

		function b = compute_uStar(self)

			% check variables KLUDGE! UPDATE WHEN YOU MAKE THE VECTORS TIME VARYING OR NOT
			uStar     = self.coefficients.uStar;
			r         = self.coefficients.r;
			if Coefficients.checkFunctionForVariable(uStar,'t') == 0
				uStar = @(x1,x2,t)(self.coefficients.uStar(x1,x2));
			end
			if Coefficients.checkFunctionForVariable(r,'t') == 0
				r = @(x1,x2,t)(self.coefficients.r(x1,x2));
			end
			%---------------------------------------------------------------------------%

			% unpack variables
			nNodes    = size(self.domain.Mesh.Nodes,2);
			nElem3    = size(self.domain.Mesh.Elements,2);
			coords    = self.domain.Mesh.Nodes';
			elements3 = self.domain.Mesh.Elements';
			%uStar     = self.coefficients.uStar;
			%r         = self.coefficients.r;

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

	end
end
