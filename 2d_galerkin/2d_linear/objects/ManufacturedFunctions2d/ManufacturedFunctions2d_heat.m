classdef ManufacturedFunctions2d_heat < ManufacturedFunctions2d_parabolic
% ManufacturedFunctions2d_heat(p,k,uTrue) manufactures the source for a heat equation MMS test.
%
% The heat equation is 
%
%			(pu)_t - div (k grad u) = f
%
% The ManufacturedFunctions2d_heat subclass manufactures the source term f for
% this problem. 
%
% Author: Tyler Fara			Date: July 28, 2024
%-----------------------------------------------------------------------------%
% Notes
%  (1) The heat equation is the simplest parabolic PDE; all parabolic PDEs will
%  have a p and k coefficient. Thus, the ManufacturedFunctions2d superclass
%  stores coefficients p and k as properties. Since the heat equation has no
%  additional coefficients, the ManufacturedFunctions2d_heat subclass has no
%  additional properties (other than the properties already stored by the
%  ManufacturedFunctions2d superclass).
%-----------------------------------------------------------------------------%


	properties
	end

	methods
		function self = ManufacturedFunctions2d_heat(c,k,uTrue,NameValueArgs)
		% ManufacturedFunctions2d_heat(c,k,uTrue) inputs are symfun objects
			
			arguments
				c
				k
				uTrue
				NameValueArgs.u_N = 0;
				NameValueArgs.alpha_R = 0;
				NameValueArgs.u_R = 0;
			end
			
			% call superclass constructor
			self@ManufacturedFunctions2d_parabolic(c,k,uTrue,NameValueArgs);

			% manufacture RHS
			self.f = self.manufactureRHS;

		end

		function f = manufactureRHS(self)
		% Manufactures symfun f (i.e. the right-hand side) from data

			% manufacture RHS
			x = sym('x',[1 2]); syms t;
			f = self.cu_t + self.divq;

		end

		function funcs = functionHandles(self)

			% call superclass method
			funcs = functionHandles@ManufacturedFunctions2d_parabolic(self);

		end

	end

end
