classdef ManufacturedFunctions2d_heat < ManufacturedFunctions2d
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
		function self = ManufacturedFunctions2d_heat(p,k,uTrue)
		% ManufacturedFunctions2d_heat(p,k,uTrue) inputs are symfun objects
			
			% call superclass constructor
			self@ManufacturedFunctions2d(p,k,uTrue)

		end

		function f = manufactureRHS(self)
		% Manufactures symfun f (i.e. the right-hand side) from data

			% manufacture RHS
			x = sym('x',[1 2]); syms t;
			f = self.p * self.u_t + self.divq;

		end


		function cofs = coefficients2FunctionHandles(self)
		% Converts symfun coefficients p and k to structure containing function_handles
		%	For cofs = self.coefficients2FunctionHandles outputs will be of
		%	format cofs.k and cofs.r

			% convert coefficients to function handles, store as struct
			cofs.k = matlabFunction(self.k);
			cofs.p = matlabFunction(self.p);

		end

	end

end
