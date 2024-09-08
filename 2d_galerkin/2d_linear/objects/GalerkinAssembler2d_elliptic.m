classdef GalerkinAssembler2d_elliptic < GalerkinAssembler2d

	properties
	end

	methods
		function self = GalerkinAssembler2d_elliptic
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
			x = sym('x',[1 2],'real');
			vars = x;

			% if inputs passed by boundary condition type
			if nargin >= 6

				% check function variables
				u_D = GalerkinAssembler2d_parabolic.getFunctionHandles(varargin{1},vars);
A				u_N = GalerkinAssembler2d_parabolic.getFunctionHandles(varargin{2},vars);
				u_R = GalerkinAssembler2d_parabolic.getFunctionHandles(varargin{3},vars);
				alpha = GalerkinAssembler2d_parabolic.getFunctionHandles(varargin{4},vars);

				% call parent method
				dom = GalerkinAssembler2d.assembleBoundary(dom,bcTypes,0,u_N,0,0);
			
			% if inputs passed by edge
			elseif nargin == 3

				% check function variables
				x = sym('x',[1 2],'real');;
				vars = x;
				bcConds = GalerkinAssembler2d_parabolic.getFunctionHandles(varargin{1},vars);

				% call parent method	
				dom = GalerkinAssembler2d.assembleBoundary(dom,bcTypes,bcConds);
			end
		end
	end
end
