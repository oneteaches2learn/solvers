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
			
			elseif nargin == 11
				bcTypes_inc = varargin{9};
				bcTypes = [bcTypes bcTypes_inc];
			end

			% create symbolic variables
			x = sym('x',[1 2],'real'); syms t;
			vars = [x t];

			% if inputs passed by boundary condition type
			if nargin == 6 || nargin == 7

				% check function variables
				u_D = GalerkinAssembler2d_parabolic.getFunctionHandles(varargin{1},vars);
				u_N = GalerkinAssembler2d_parabolic.getFunctionHandles(varargin{2},vars);
				alpha = GalerkinAssembler2d_parabolic.getFunctionHandles(varargin{3},vars);
				u_R = GalerkinAssembler2d_parabolic.getFunctionHandles(varargin{4},vars);

				% call parent method
				dom = GalerkinAssembler2d.assembleBoundary(dom,bcTypes,u_D,u_N,alpha,u_R);
			
			% if inputs passed by boundary condition type, including dynamic boundary conditions
			elseif nargin == 10 || nargin == 11

				% check function variables
				u_D = GalerkinAssembler2d_parabolic.getFunctionHandles(varargin{1},vars);
				u_N = GalerkinAssembler2d_parabolic.getFunctionHandles(varargin{2},vars);
				alpha_R = GalerkinAssembler2d_parabolic.getFunctionHandles(varargin{3},vars);
				u_R = GalerkinAssembler2d_parabolic.getFunctionHandles(varargin{4},vars);
				alpha_T = GalerkinAssembler2d_parabolic.getFunctionHandles(varargin{5},vars);
				beta_T = GalerkinAssembler2d_parabolic.getFunctionHandles(varargin{6},vars);
				gamma_T = GalerkinAssembler2d_parabolic.getFunctionHandles(varargin{7},vars);
				u_T = GalerkinAssembler2d_parabolic.getFunctionHandles(varargin{8},vars);

				% call parent method
				dom = GalerkinAssembler2d.assembleBoundary( ...
							dom,bcTypes,u_D,u_N,alpha_R,u_R,alpha_T,beta_T,gamma_T,u_T);

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

		function auxfun = assembleCoefficients(c,k,r,f,uInit)

			% call parent method
			auxfun = GalerkinAssembler2d.assembleCoefficients(k,r,f);

			% store time-varying coefficient
			x = sym('x',[1 2],'real'); syms t;
			auxfun.cofs.c = matlabFunction(symfun(c,x));

			% store initial condition
			auxfun.uInit = matlabFunction(symfun(uInit,x));

		end

		function dom = assembleTimeStepping(dom,T,dt,eq)

			% create time-stepping object
			if nargin == 3
				time = TimeStepping(T,dt);
			elseif nargin == 4
				time = TimeStepping(T,dt,eq);
			end

			% set time-stepping mesh
			time = time.setMesh;

			% store result
			dom.time = time;

		end

	end
end
