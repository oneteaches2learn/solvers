classdef GalerkinAssembler2d_parabolic < GalerkinAssembler2d

	properties
	end

	methods
		function self = GalerkinAssembler2d_parabolic
		end
	end

	methods (Static)
		function dom = assembleBoundary(dom,bcTypes,varargin); 

			% old code split bcTypes and bcTypes_inc, but now they are combined
			if nargin == 7
				bcTypes_inc = varargin{5};
				bcTypes = [bcTypes bcTypes_inc];
			end

			% create symbolic variables
			x = sym('x',[1 2],'real'); syms t;
			vars = [x t];

			% if inputs passed by boundary condition type
			if nargin >= 6

				% check function variables
				u_D = GalerkinAssembler2d_parabolic.getFunctionHandles(varargin{1},vars);
				u_N = GalerkinAssembler2d_parabolic.getFunctionHandles(varargin{2},vars);
				u_R = GalerkinAssembler2d_parabolic.getFunctionHandles(varargin{3},vars);
				alpha = GalerkinAssembler2d_parabolic.getFunctionHandles(varargin{4},vars);

				% call parent method
				dom = GalerkinAssembler2d.assembleBoundary(dom,bcTypes,u_D,u_N,alpha,u_R);
			
			% if inputs passed by edge
			elseif nargin == 3

				% check function variables
				x = sym('x',[1 2],'real'); syms t;
				vars = [x t];
				bcConds = GalerkinAssembler2d_parabolic.getFunctionHandles(varargin{1},vars);

				% call parent method	
				dom = GalerkinAssembler2d.assembleBoundary(dom,bcTypes,bcConds);
			end
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
