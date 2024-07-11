classdef GalerkinHeat2d_mms
% DOCUMENTATION NEEDED! In the meantime, see the documentation for
% poissonFD1d_finegrid, which is a very similar object. 
%
% Notes
%	(1) Use the member funtion self.instantiateFineGridTest to create a fine
%	grid test using the same input data as the corresponding mms test. This can
%	be used to compare results for an mms versus fine grid test.
%-----------------------------------------------------------------------------%

	properties
		problems
		auxFunctions
		domain
		timeStepping
		errType
		errors
		ratios
		orders
	end

	methods
		function self = GalerkinHeat2d_mms(dom,time,auxfun,NameValueArgs)
			arguments
				dom 	Domain2d
				time	TimeStepping
				auxfun  HeatAuxFunctions2d
				NameValueArgs.errType = "L2"
				NameValueArgs.demo = 0
			end

			% store inputs
			self.auxFunctions = auxfun;
			self.errType = NameValueArgs.errType;
			self.timeStepping = time;

			% manufacture boundary conditions
			fprintf(' Setting BCs:'), tic
			self.domain = self.manufactureBoundaryConditions(dom,auxfun);
			executionTime = toc;
			fprintf(' %f s\n',executionTime)

			if NameValueArgs.demo == 0
				self.problems = self.solveManufacturedProblems(NameValueArgs.demo);
				[self.errors,self.ratios,self.orders] = self.computeErrors;
			else
				self.problems = self.solveManufacturedProblems(NameValueArgs.demo);
			end

		end


		function problems = solveManufacturedProblems(self,demo)

			% Convert to function_handles
			f = self.auxFunctions.source2FunctionHandle;
			cofs = self.auxFunctions.coefficients2FunctionHandles;
			uInit = self.auxFunctions.uInit2FunctionHandle;

			% Solve problems for successive base^-p
			fprintf('MMS Test Begun\n')
			fprintf('Solving Problems\n')

			base = 2; % <~~ at some point, need to pass this instead of hard-setting it

			% run mms test
			if demo == 0, p_min = 1; p_max = 4; else p_min = demo; p_max = demo; end
			for p = p_min:p_max

				tic 
				% successively refine mesh and time-stepping 
				dom_p  = self.domain.setMesh(p,base);
				time_p = self.timeStepping.setMesh(2*(p-1),base); %<~~ NOTE: time stepping one order higher

				% run solver on current mesh
				prob_p = GalerkinHeat2d_solver(dom_p,time_p,cofs,uInit,f);

				% store results
				if demo == 0, problems{p} = prob_p;
				else problems{1} = prob_p;
				end

				executionTime = toc;
				fprintf(' p = %i solved: %f s\n',p,executionTime)
			end
		end

		function dom = manufactureBoundaryConditions(self,dom,auxfun);

			% unpack coefficients
			x = sym('x',[1 2]); syms t;
			uTrue = self.auxFunctions.uTrue;
			q = self.auxFunctions.q;
			nEdges = dom.NumEdges;

			% Manufacture Dirichlet BC
			u_d = uTrue;
			u_d = matlabFunction(u_d);

			% assign boundary functions
			for i = 1:nEdges

				% assign Dirichlet BC
				if dom.edges(i).boundaryType == 'D'
					dom.edges(i).boundaryCondition = u_d;

				% assign Neumann BC
				elseif dom.edges(i).boundaryType == 'N'
					n_i = dom.edges(i).outwardNormal;
					g_i = symfun(dot(q,n_i),[x t])
					g_i = matlabFunction(g_i);
					dom.edges(i).boundaryCondition = g_i;

				% assign Robin BC
				elseif dom.edges(i).boundaryType == 'R'
					alpha_i = symfun(1.0,[x t]);
					n_i = dom.edges(i).outwardNormal;
					g_i = symfun(uTrue - dot(q,n_i) / alpha_i,[x t]);
					alpha_i = matlabFunction(alpha_i);
					g_i = matlabFunction(g_i);
					dom.edges(i).boundaryCondition = {alpha_i,g_i};
				end
			end

		end


		function [errors,ratios,orders] = computeErrors(self)

			% store variables
			trials = length(self.problems);
			base = self.problems{1}.domain.base;

			% Store errors
			errors = zeros(1,trials);
			ratios = zeros(1,trials-1);
			orders = zeros(1,trials-1);

			% Compute error
			errorObj = mmsErrorComputer(self);
			errors = errorObj.errors;

			% Compute ratios and orders
			for i = 2:trials
				ratios(i-1) = errors(i-1)/errors(i);
				orders(i-1) = log(ratios(i-1)) / log(base);
			end
		end


		function v = computeOutwardNormal(self,nod1,nod2)
			x = nod2(1) - nod1(1);
			y = nod2(2) - nod1(2);
			v = [y,-x];
		end


		function ax = plot(self)

			h = zeros(1,4);
			for i = 1:4
				base = self.problems{i}.domain.base;
				p = self.problems{i}.domain.p;
				h(i) = base^-p;
			end

			loglog(h,h,'k--',h,h.^2,'k-.',h,self.errors,'r*-','LineWidth',3,'MarkerSize',18); grid on;
			title('Error of D_hf');
			legend('linear','quadratic','e(h)','location','southeast');
			f = gcf;
			f.Position = [100 100 800 * 3/4 600 * 3/4];

			ax = gca;

		end

		function fineGrid = instantiateFineGridObject(self)
	
			dom = self.problems{1}.domain;
			bc  = self.problems{1}.boundary;
			k   = self.problems{1}.conductivity;
			f   = self.problems{1}.source;

			fineGrid = poissonFD1d_finegrid(dom,bc,k,f);

		end

		function latexTable(self,fid)
			
			aType = self.problems{1}.boundary.aType;
			bType = self.problems{1}.boundary.bType;
			errors = self.errors;
			orders = self.orders;

			% Table Preamble
			fprintf(fid,'%%\n')
			fprintf(fid,'\\begin{table}[]\n')
			fprintf(fid,'\t\\centering\n')
			fprintf(fid,'\t\\begin{tabular}{c c c}\n')

			% Table title line
			fprintf(fid,'\t\\hline\n');
			fprintf(fid,'\t\t refinement & error & rate \\\\ \n');
			fprintf(fid,'\t\\hline\n');

			% Table Main Data
			fprintf(fid,'\t\t $h$ & %.4d & --     \\\\ \n',errors(1));
			for i = 1:3
				fprintf(fid,'\t\t $h^%i$ & %.4d & %.4f \\\\ \n',i+1,errors(i+1),orders(i));
			end
			fprintf(fid,'\t\\hline\n')

			% Table Footer
			fprintf(fid,'\t\\end{tabular}\n')
			fprintf(fid,'\t\\caption{\\textbf{Problem :} MMS, Boundary: %c and %c}\n',aType,bType)
			fprintf(fid,'\t\\label{tab:1-}\n')
			fprintf(fid,'\\end{table}\n')
			fprintf(fid,'%%\n')
		end

		function L2_IP = L2_IP(self,prob,timestep,quadOrder)

			t = timestep * prob.time.dt;

			% set quadrature order
			if nargin < 4, quadOrder = 4; end

			% store variables
			nElem = size(prob.domain.Mesh.Elements,2);
			uTrue = matlabFunction(self.uTrue);
			uTrue = @(x1,x2)(uTrue(x1,x2,t));
			u_h   = prob.solution(:,timestep);

			% loop over elements
			L2_IP = 0;
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

		function err = L2_L2_error(self,prob,quadOrder)

			if nargin < 3, quadOrder = 4; end
			
			dt  = prob.time.dt;
			N_t = prob.time.N_t;

			err = 0;
			for i = 1:N_t+1
				err = err + dt * self.L2_IP(prob,i,quadOrder);
			end

			err = sqrt(err);
		end

		function err = Linfty_L2_error(self,prob,quadOrder)

			if nargin < 3, quadOrder = 4; end
			
			dt  = prob.time.dt;
			N_t = prob.time.N_t;

			err = 0;
			for i = 1:N_t+1
				err_i = sqrt(self.L2_IP(prob,i,quadOrder));
				if err_i > err
					err = err_i;
				end
			end

		end


		function plotTrueSolution(self)

			nProb = length(self.problems);
			prob = self.problems{nProb};
			uTrue = self.uTrue;
			uTrue = matlabFunction(uTrue);

			% plot final timestep unless otherwise specified
			if nargin < 2, timestep = prob.time.N_t; end
			t = timestep * prob.time.dt;

			% store domain information
			coordinates = prob.domain.Mesh.Nodes';
			elements3 = prob.domain.Mesh.Elements';
			elements4 = [];

			% get solution at final time step
			U = uTrue(coordinates(:,1),coordinates(:,2),t * ones(size(coordinates,1),1))

			% plot data
			trisurf(elements3,coordinates(:,1),coordinates(:,2),U', ...
				'facecolor','interp')
			hold on
			trisurf(elements4,coordinates(:,1),coordinates(:,2),U', ...
				'facecolor','interp')
			hold off

			% format plot
			view(10,40);
			title('Solution of the Problem')
		end


	end
end

