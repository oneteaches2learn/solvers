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
				alpha = GalerkinAssembler2d_parabolic.getFunctionHandles(varargin{3},vars);
				u_R = GalerkinAssembler2d_parabolic.getFunctionHandles(varargin{4},vars);

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
