classdef GalerkinParabolic2dSemilinear_solver < GalerkinParabolic2d_solver

	properties
	end

	properties (Hidden)
		iter 
		tolerance = 1e-6
	end

	properties (Hidden,Dependent)
		isFirstIter
		stoppingCriterion
	end

	methods
		function self = GalerkinParabolic2dSemilinear_solver(dom,time,cofs,uInit,f)
			
			self@GalerkinParabolic2d_solver(dom,time,cofs,uInit,f);

		end

		function self = solve(self)

			% initialize problem
			FreeNodes = self.domain.freeNodes;
			self = self.initializeProblem;

			% loop on timesteps
			for timestep = 1:self.time.M_t

				self.timestep = timestep;
				self = self.initializeTimestep;

				% resolve potential nonlinearity
				while self.stoppingCriterion == 0

					% prepare iteration
					self = self.initializeIteration;
					
					% assemble problem
					self  = self.assembleTensors;
					self  = self.assembleVectors;
					self  = self.assembleBCs;
					[S,b] = self.finalAssembly;

					% solve and store solution
					v = sparse(self.domain.nNodes,1);
					v(FreeNodes) = S(FreeNodes,FreeNodes) \ b(FreeNodes);
					self.vectors.U_currIter = v + self.vectors.U_D;

				end

				% store result
				self.solution(:,self.timestep) = self.vectors.U_currIter;

			end

			self = self.cleanup;

		end


		% INITIALIZATION FUNCTIONS -------------------------------------------%
		function self = initializeVectors(self)

			% call superclass method
			self = initializeVectors@GalerkinParabolic2d_solver(self);

			% create fields for vector storage
			self.vectors.U_prevIter = [];
			self.vectors.U_currIter = [];

		end

		function self = initializeTimestep(self)

			% call superclass method
			self = initializeTimestep@GalerkinParabolic2d_solver(self);

			% update vectors for next iteration
			self.vectors.U_prevIter = self.vectors.U_prevTime;
			self.vectors.U_currIter = [];

			% reset iteration counter
			self.iter = 0;

		end

		function self = initializeIteration(self)

			% update iteration counter
			self.iter = self.iter + 1;

			% if past first iteration, shuffle vectors for next iteration
			if self.iter ~= 1
				self.vectors.U_prevIter = self.vectors.U_currIter;
				self.vectors.U_currIter = [];
			end

		end


		% ASSEMBLY FUNCTIONS -------------------------------------------------%
		function self = assembleTensors(self)

			% if first timestep, create tensors
			if self.isFirstTimeStep == 1 && self.isFirstIter == 1
				self.tensors.A   = self.assembleStiffnessMatrix;
				self.tensors.M_p = self.assembleMassMatrix(self.coefficients.p);

			% else, update tensors as needed
			else

				% update A
				if self.tensors.timeVarying.A == 1 && self.isFirstIter == 1
					self.tensors.A = self.assembleStiffnessMatrix;
				end

				% update M_p
				if self.tensors.timeVarying.M_p == 1 && self.isFirstIter == 1
					self.tensors.M_p_prevTime = self.tensors.M_p;
					cof = self.coefficients.p;
					self.tensors.M_p = self.assembleMassMatrix(cof);
				end
			end

		end

		function self = assembleVectors(self)

			% if first timestep, create vectors
			if self.isFirstTimeStep == 1
				self.vectors.b_vol = self.computeVolumeForces;

			% else, update vectors as needed
			else

				% update b_vol
				if self.vectors.timeVarying.b_vol == 1
					self.vectors.b_vol = self.computeVolumeForces;
				end

			end

		end

		function A = assembleStiffnessMatrix(self)

			% store variables
			k = self.coefficients.k;

			% check coefficient variables
			if Coefficients.isTimeVarying(k) == 0
				k = @(x1,x2,t)(k(x1,x2));
			end

			% unpack variables
			nNodes    = size(self.domain.Mesh.Nodes,2);
			nElem3    = size(self.domain.Mesh.Elements,2);
			coords    = self.domain.Mesh.Nodes';
			elements3 = self.domain.Mesh.Elements';

			% initialize storage
			A = sparse(nNodes,nNodes);

			% compute k on nodes
			for j = 1:nNodes
				K(j) = k(coords(j,1),coords(j,2),self.t);
			end

			% assemble stiffness matrix
			for j = 1:nElem3
				elementInd   = elements3(j,:);
				elementCoord = coords(elementInd,:);
				K_j = (1/3) * sum(K(elementInd));
				A(elementInd,elementInd) = A(elementInd,elementInd) + ...
					K_j * self.stima3(elementCoord);
			end

		end

		function B = assembleMassMatrix(self,c)

			% check coefficient variables
			if Coefficients.isTimeVarying(c) == 0
				c = @(x1,x2,t)(c(x1,x2));
			end
				
			% unpack variables
			nNodes    = size(self.domain.Mesh.Nodes,2);
			nElem3    = size(self.domain.Mesh.Elements,2);
			coords 	  = self.domain.Mesh.Nodes';
			elements3 = self.domain.Mesh.Elements';

			% initialize storage
			B = sparse(nNodes,nNodes);
			
			% assemble base mass matrix
			for j = 1:nElem3
				elementInd = elements3(j,:);
				elementCoord = coords(elementInd,:);
				B(elementInd,elementInd) = B(elementInd,elementInd) + ...
					self.domain.elemAreas(j) * [2,1,1;1,2,1;1,1,2] / 12;
			end

			% compute r on nodes
			for j = 1:nNodes
				C(j) = c(coords(j,1),coords(j,2),self.t);
			end

			% scale mass matrix by r values
			B = B.*C';

		end

		function b = computeVolumeForces(self)

			% store variables
			f = self.f;

			% check coefficient variables
			if Coefficients.isTimeVarying(f) == 0
				f = @(x1,x2,t)(f(x1,x2));
			end

			% unpack variables
			nNodes    = size(self.domain.Mesh.Nodes,2);
			nElem3    = size(self.domain.Mesh.Elements,2);
			coords    = self.domain.Mesh.Nodes';
			elements3 = self.domain.Mesh.Elements';

			% initialize storage
			b = sparse(nNodes,1);

			% compute volume forces
			for j = 1:nElem3
				elementInd    = elements3(j,:);
				elementCoord  = coords(elementInd,:);
				b(elementInd) = b(elementInd) + ...
					self.domain.elemAreas(j) * ...
					f(sum(elementCoord(:,1))/3,sum(elementCoord(:,2))/3,self.t) / 3;
			end

		end

		% GETTERS ------------------------------------------------------------%
		function val = get.isFirstIter(self)

			val = (self.iter == 1);

		end

		function val = get.stoppingCriterion(self)

			% if entering first iteration, proceed
			if self.iter == 0
				
				val = 0;

			% else check L2 norm of difference between successive iterations
			else

				% compute L2 error
				iterDiff = self.vectors.U_currIter - self.vectors.U_prevIter;
				iterErr  = self.domain.L2norm_piecewiseLinear(iterDiff);

				% if less then tolerance, break iteration
				if iterErr < self.tolerance
					val = 1;

				% else proceed
				else
					val = 0;
				end
			end

		end

	end
end
