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

		function tensors = initializeTensors(self,t)
			
			% create fields for tensor storage
			tensors.A   = [];
			tensors.M_r = [];
			tensors.M_p = [];
			tensors.E   = [];

			% check which tensors are time-varying
			tensors.timeVarying.A   = Coefficients.checkTimeVarying(self.coefficients.k);
			tensors.timeVarying.M_r = Coefficients.checkTimeVarying(self.coefficients.r);
			tensors.timeVarying.M_p = Coefficients.checkTimeVarying(self.coefficients.p);

		end

		function self = assembleTensors(self,t)

			% if first timestep, create tensors
			if t == 1 * self.time.dt
				self.tensors.A   = self.assembleStiffnessMatrix(t);
				self.tensors.M_r = self.assembleMassMatrix(self.coefficients.r,t);
				self.tensors.M_p = self.assembleMassMatrix(self.coefficients.p,t);

			% else, update tensors as needed
			else

				% update A
				if self.tensors.timeVarying.A == 1
					self.tensors.A = self.assembleStiffnessMatrix(t);
				end

				% update M_r
				if self.tensors.timeVarying.M_r == 1
					cof = self.coefficients.r;
					self.tensors.M_r = self.assembleMassMatrix(cof,t);
				end

				% update M_p
				if self.tensors.timeVarying.M_p == 1
					cof = self.coefficients.p;
					self.tensors.M_p = self.assembleMassMatrix(cof,t);
				end

			end

			% update Robin boundary tensor
			[self.tensors.E,temp] = self.computeRobinBCs(t);

		end

		function vectors = assembleVectors(self,t)

			% assemble vectors
			vectors.uStar = self.compute_uStar(t);
			vectors.b_vol = self.computeVolumeForces(t);
			vectors.U_D   = self.computeDirichletBCs(t);
			vectors.b_neu = self.computeNeumannBCs(t);
			[temp,vectors.b_rob] = self.computeRobinBCs(t);

		end

		function [S,b] = finalAssembly(self,vectors,U_prev)

			% store variables
			tensors = self.tensors;

			% assemble LHS
			S = self.time.dt * (tensors.A + tensors.M_r + tensors.E) + tensors.M_p;

			% assemble RHS
			b = self.time.dt * (vectors.b_vol - vectors.b_neu + vectors.b_rob + ...
					vectors.uStar) - S * vectors.U_D + tensors.M_p * U_prev;

		end

		%{
		function vec = compute_uStar(self,t)

			% store variables
			nNodes = self.domain.nNodes;

			% loop over nodes
			vec = zeros(nNodes,1);
			for i = 1:nNodes
				vec(i) = self.coefficients.uStar(vec(i),vec(i),t);
			end

		end
		%}

		function b = compute_uStar(self,t)


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
					r(sum(elementCoord(:,1))/3,sum(elementCoord(:,2))/3,t) * ...
					uStar(sum(elementCoord(:,1))/3,sum(elementCoord(:,2))/3,t) / 6;

			end

		end

	end
end
