classdef GalerkinParabolic2d_assembler

	properties
	end

	methods
		function self = GalerkinParabolic2d_assembler
		end
	end

	methods (Static)
		function dom = assembleDomainGeometry(xLim,yLim,inc,eps)

			% create domain object
			if nargin == 2
				dom = Domain2d(xLim,yLim);
			elseif nargin == 4
				dom = Domain2d_punctured(xLim,yLim,inc,eps);
			end

		end

		function dom = assembleBoundary(dom,bcTypes,u_D,u_N,alpha,u_R,bcTypes_inc); 

			% check function variables
			x = sym('x',[1 2],'real'); syms t;
			u_D = GalerkinParabolic2d_assembler.getFunctionHandles(u_D);
			u_N = GalerkinParabolic2d_assembler.getFunctionHandles(u_N);
			u_R = GalerkinParabolic2d_assembler.getFunctionHandles(u_R);
			alpha = GalerkinParabolic2d_assembler.getFunctionHandles(alpha);

			% set boundary conditions on outer boundary
			for i = 1:4
				if bcTypes{i} == 'D'
					if iscell(u_D)
						bcConds{i} = u_D{i};
					else
						bcConds{i} = u_D;
					end
				elseif bcTypes{i} == 'N'
					if iscell(u_N)
						bcConds{i} = u_N{i};
					else
						bcConds{i} = u_N;
					end
				elseif bcTypes{i} == 'R'
					if iscell(alpha) && ~iscell(u_R)
						bcConds{i} = {alpha{i},u_R};
					elseif ~iscell(alpha) && iscell(u_R)
						bcConds{i} = {alpha,u_R{i}};
					elseif iscell(alpha) && iscell(u_R)
						bcConds{i} = {alpha{i},u_R{i}};
					else
						bcConds{i} = {alpha,u_R};
					end
				end
			end

			% set boundary conditions on inclusions
			if nargin == 7
				if bcTypes_inc == 'D'
					if iscell(u_D)
						bcConds_inc = u_D{5};
					else
						bcConds_inc = u_D;
					end
				elseif bcTypes_inc == 'N'
					if iscell(u_N)
						bcConds_inc = u_N{5};
					else
						bcConds_inc = u_N;
					end
				elseif bcTypes_inc == 'R'
					if iscell(alpha) && ~iscell(u_R)
						bcConds_inc = {alpha{5},u_R};
					elseif ~iscell(alpha) && iscell(u_R)
						bcConds_inc = {alpha,u_R{5}};
					elseif iscell(alpha) && iscell(u_R)
						bcConds_inc = {alpha{5},u_R{5}};
					else
						bcConds_inc = {alpha,u_R};
					end
				end
			end

			% create Boundary2d object
			if nargin == 7
				bound = Boundary2d_punctured(bcTypes,bcConds,bcTypes_inc,bcConds_inc);
			else
				bound = Boundary2d(bcTypes,bcConds);
			end

			% set boundary conditions on domain
			dom = dom.setEdgeBCTypes(bound);
			dom = dom.setEdgeBCConditions(bound);

		end

		function f = getFunctionHandles(f)

			% check function variables
			x = sym('x',[1 2],'real'); syms t;
			if iscell(f)
				for i = 1:length(f)
					f{i} = symfun(f{i},[x t]);
				end
			else
				f = symfun(f,[x t]);
			end

			% convert to function handles
			if iscell(f)
				for i = 1:length(f)
					f{i} = matlabFunction(f{i});
				end
			else
				f = matlabFunction(f);
			end

		end

		function dom = assembleMesh(dom,p,base)

			% set mesh on domain
			dom = dom.setMesh(p,base);

		end

		function cofs = assembleCoefficients(p,k)

			% check function variables
			x = sym('x',[1 2],'real'); syms t;
			p = symfun(p,x);
			k = symfun(k,x);

			% convert to function_handles
			cofs.p = matlabFunction(p);
			cofs.k = matlabFunction(k);

		end

		function uInit = assembleInitialCondition(uInit)

			% check function variables
			x = sym('x',[1 2],'real'); syms t;
			uInit = symfun(uInit,[x t]);

			% convert to function_handles
			uInit = matlabFunction(uInit);

		end

		function time = assembleTimeStepping(T,dt,eq)

			% create time-stepping object
			if nargin == 2
				time = TimeStepping(T,dt);
			elseif nargin == 3
				time = TimeStepping(T,dt,eq);
			end

			% set time-stepping mesh
			time = time.setMesh;

		end

		function f = assembleSource(f)

			% check function variables
			x = sym('x',[1 2],'real'); syms t;
			f = symfun(f,x);

			% convert to function_handle
			f = matlabFunction(f);

		end

	end

	

end
