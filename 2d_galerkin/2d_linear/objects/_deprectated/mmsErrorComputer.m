classdef mmsErrorComputer

	properties
		mmsObj
		errors
	end

	methods
		function self = mmsErrorComputer(mmsObj)

			self.mmsObj = mmsObj;
			pmin = mmsObj.mmsParams.pmin;

			% determine number of trials stored in mmsObj
			trials = length(mmsObj.problems);

			% compute and store errors
			fprintf('Computing Errors\n');
			self.errors = zeros(1,trials);
			for i = 1:trials
				fprintf(' p = %i computed: ',i + pmin - 1);
				self.errors(i) = self.computeError(i);
			end

		end

		function err = computeError(self,trial_num)

			% get appropriate problem
			prob = self.mmsObj.problems{trial_num};
			errType = self.mmsObj.errType;
			quadOrder = 1; % <~~quad order hard coded
			
			% send problem to appropriate error module
			tic 
			% check for stationary error type
			if strcmp(errType,"L2")
				err = self.L2_error(prob,quadOrder);

			elseif strcmp(errType,"Linfty")

			elseif strcmp(errType,"L2(L2)")
				err = self.L2_L2_error(prob,quadOrder);

			elseif strcmp(errType,"Linfty(L2)")
				err = self.Linfty_L2_error(prob,quadOrder);

			else 
				error("Invalid errType")

			end

			executionTime = toc;
			fprintf('%f s\n',toc);

		end


		function err = L2_error(self,prob,quadOrder)

			% compute error
			err = 0;
			err = self.L2_IP(prob,quadOrder);
			err = sqrt(err);

		end


		function err = Linfty_error(self,prob,quadOrder)

			% compute error
			err = 0;
			err = max(prob.solution);

		end


		function err = L2_L2_error(self,prob,quadOrder)

			% store variables
			dt  = prob.time.dt;
			N_t = prob.time.N_t;

			% compute error
			err = 0;
			wait = waitbar(0,'Computing Error, Please Wait');
			for i = 1:N_t
				err = err + dt * self.L2_IP(prob,quadOrder,i);
				waitbar(i/(N_t+1),wait,'Computing Error, Please Wait');
			end
			err = sqrt(err);
			close(wait);
		end

		function err = Linfty_L2_error(self,prob,quadOrder)

			% store variables	
			dt  = prob.time.dt;
			N_t = prob.time.N_t;

			% compute error
			err = zeros(1,N_t);

			% for small N_t, parallel processing adds time
			if N_t < 100
				for i = 1:N_t
					err(i) = sqrt(self.L2_IP(prob,quadOrder,i));
				end

			% else, use parallel processing
			else
				parfor i = 1:N_t
					err(i) = sqrt(self.L2_IP(prob,quadOrder,i));
				end

			end

			err = max(err);

		end

		function L2_IP = L2_IP(self,prob,quadOrder,timestep)

			if quadOrder == 1
				L2_IP = L2_IP_firstOrder(self,prob,quadOrder,timestep);
			else
				L2_IP = L2_IP_higherOrder(self,prob,quadOrder,timestep);
			end

		end

		function L2_IP = L2_IP_firstOrder(self,prob,quadOrder,timestep)

			% check if timestep passed
			if nargin == 4, timeVarying = 1;
			else timeVarying = 0;
			end

			% if timetep passed, set current time
			if timeVarying == 1, t = (timestep-1) * prob.time.dt; end

			% store uTrue at current time
			uTrue = matlabFunction(self.mmsObj.auxFunctions.uTrue);
			if timeVarying == 1, uTrue = @(x1,x2)(uTrue(x1,x2,t)); end
			UTrue = uTrue(prob.domain.mesh.nodes(:,1),prob.domain.mesh.nodes(:,2));

			if sum(size(UTrue)) == 2
				UTrue = UTrue * ones(size(prob.solution,1),1);
			end

			% store numerical solution
			if timeVarying == 1, u_h = prob.solution(:,timestep);
			else u_h = full(prob.solution); end

			L2_IP = prob.domain.L2_IP_nodalQuadrature(UTrue-u_h,UTrue-u_h);

			%{
			% loop over elements
			L2_IP = 0;
			nElem = prob.domain.nElem; 
			for i = 1:nElem

				% interpolate u_h locally 
				locNodes = prob.domain.Mesh.Elements(:,i);

				% compute local error function on quad points
				uTrue_Qp = UTrue(locNodes);
				uSol_Qp  = u_h(locNodes); 
				err_Qp   = uTrue_Qp - uSol_Qp; 

				% quadrature on local error function
				L2_IP = L2_IP + prob.domain.elemAreas(i) * sum(err_Qp.^2) / 3; 

			end
			%}

		end

		function L2_IP = L2_IP_higherOrder(self,prob,quadOrder,timestep)

			% check if timestep passed
			if nargin == 4, timeVarying = 1;
			else timeVarying = 0;
			end

			% if timetep passed, set current time
			if timeVarying == 1, t = (timestep-1) * prob.time.dt; end

			% store uTrue at current time
			uTrue = matlabFunction(self.mmsObj.auxFunctions.uTrue);
			if timeVarying == 1, uTrue = @(x1,x2)(uTrue(x1,x2,t)); end

			% store numerical solution
			if timeVarying == 1, u_h   = prob.solution(:,timestep);
			else u_h = full(prob.solution); end

			% loop over elements
			L2_IP = 0;
			nElem = size(prob.domain.mesh.Mesh.Elements,2);
			for i = 1:nElem
				% interpolate u_h locally 
				locNodes = prob.domain.mesh.Mesh.Elements(:,i); 
				locCoord = prob.domain.mesh.Mesh.Nodes(:,locNodes)';
				uLoc = scatteredInterpolant(locCoord,u_h(locNodes)); 

				% compute local error function on quad points
				Qp = quadtriangle(quadOrder,'Domain',locCoord); 
				uTrue_Qp = uTrue(Qp.Points(:,1),Qp.Points(:,2)); 
				uSol_Qp  = uLoc(Qp.Points(:,1),Qp.Points(:,2)); 
				err_Qp   = uTrue_Qp - uSol_Qp; 

				% quadrature on local error function
				L2_IP = L2_IP + dot(Qp.Weights, err_Qp .* err_Qp); 
			end
		end

	end

end

