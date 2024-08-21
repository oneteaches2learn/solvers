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
				dom = PuncturedDomain2d(xLim,yLim,inc,eps);
			end

		end

		function dom = assembleBoundary(dom,bcTypes,u_D,u_N,alpha,u_R,bcTypes_inc); 

			% check function variables
			x = sym('x',[1 2],'real'); syms t;
			u_D = symfun(u_D,[x t]);
			u_N = symfun(u_N,[x t]);
			u_R = symfun(u_R,[x t]);
			alpha = symfun(alpha,[x t]);

			% convert to function_handles
			u_D = matlabFunction(u_D);
			u_N = matlabFunction(u_N);
			u_R = matlabFunction(u_R);
			alpha = matlabFunction(alpha);

			% set boundary conditions on outer boundary
			for i = 1:4
				if bcTypes{i} == 'D'
					bcConds{i} = u_D;
				elseif bcTypes{i} == 'N'
					bcConds{i} = u_N;
				elseif bcTypes{i} == 'R'
					bcConds{i} = {alpha,u_R};
				end
			end

			% set boundary conditions on inclusions
			if nargin == 7
				if bcTypes_inc == 'D'
					bcConds_inc = u_D;
				elseif bcTypes_inc == 'N'
					bcConds_inc = u_N;
				elseif bcTypes_inc == 'R'
					bcConds_inc = {alpha,u_R};
				end
			end

			% create Boundary2d object
			if nargin == 7
				bound = PuncturedBoundary2d(bcTypes,bcConds,bcTypes_inc,bcConds_inc);
			else
				bound = Boundary2d(bcTypes,bcConds);
			end

			% set boundary conditions on domain
			dom = dom.setEdgeBCTypes(bound);
			dom = dom.setEdgeBCConditions(bound);

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
