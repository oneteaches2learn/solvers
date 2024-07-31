classdef mmsErrorComputer

	properties
		mmsObj
		errors
	end

	methods
		function self = mmsErrorComputer(mmsObj)

			self.mmsObj = mmsObj;

			% determine number of trials stored in mmsObj
			trials = length(mmsObj.problems);

			% compute and store errors
			fprintf('Computing Errors\n');
			self.errors = zeros(1,trials);
			for i = 1:trials
				fprintf(' p = %i computed: ',i);
				self.errors(i) = self.computeError(i);
			end

		end

		function err = computeError(self,trial_num)

			% get appropriate problem
			prob = self.mmsObj.problems{trial_num};
			errType = self.mmsObj.errType;
			quadOrder = 2;
			
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
			for i = 1:N_t+1
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
			wait = waitbar(0,'Computing Error, Please Wait');
			err = 0;
			for i = 1:N_t+1
				err_i = sqrt(self.L2_IP(prob,quadOrder,i));
				if err_i > err
					err = err_i;
				end
				waitbar(i/(N_t+1),wait,'Computing Error, Please Wait');
			end
			close(wait)

		end

		function L2_IP = L2_IP(self,prob,quadOrder,timestep)

			% check if timestep passed
			if nargin == 4, timeVarying = 1;
			else timeVarying = 0;
			end

			% if timetep passed, set current time
			if timeVarying == 1, t = timestep * prob.time.dt; end

			% store uTrue
			uTrue = matlabFunction(self.mmsObj.uTrue);
			if timeVarying == 1, uTrue = @(x1,x2)(uTrue(x1,x2,t)); end

			% store numerical solution
			if timeVarying == 1, u_h   = prob.solution(:,timestep);
			else u_h = full(prob.solution); end

			% loop over elements
			L2_IP = 0;
			nElem = size(prob.domain.Mesh.Elements,2);
			for i = 1:nElem
				% interpolate u_h locally 
				locNodes = prob.domain.Mesh.Elements(:,i); 
				locCoord = prob.domain.Mesh.Nodes(:,locNodes)';
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

