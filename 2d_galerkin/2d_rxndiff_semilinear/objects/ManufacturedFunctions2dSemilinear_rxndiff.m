classdef ManufacturedFunctions2dSemilinear_rxndiff < ManufacturedFunctions2d_rxndiff
% ManufacturedFunctions2dSemilinear_rxndiff(p,k,r,uTrue) manufactures the source for a semilinear reaction-diffusion equation MMS test.
%
% The semilinear reaction-diffusion equation is 
%
%			(pu)_t - div (k grad u) + r(u) = f
%
% The ManufacturedFunctions2Semilinear_rxndiff subclass manufactures the source
% term f for this problem. 
%
% Author: Tyler Fara			Date: July 31, 2024
%-----------------------------------------------------------------------------%

	properties
	end

	methods
		function self = ManufacturedFunctions2dSemilinear_rxndiff(p,k,r,uTrue)
		% ManufacturedFunctions2d_pennes(p,k,r,uTrue) inputs are symfun objects
			
			% call superclass constructor
			self@ManufacturedFunctions2d_rxndiff(p,k,r,uTrue);

			% manufacture RHS
			self.f = self.manufactureRHS; 

		end

		function f = manufactureRHS(self)
		% Manufactures symfun f (i.e. the right-hand side) from data

			% manufacture RHS
			x = sym('x',[1 2]); syms t;
			f = self.pu_t + self.divq + self.r * self.uTrue;

		end

		function cofs = coefficients2FunctionHandles(self)
		% Converts symfun coefficients p, k, r to structure containing function_handles
		%	For cofs = self.coefficients2FunctionHandles outputs will be of
		%	format cofs.k, cofs.r, cofs.p;

			% convert coefficients to function handles, store as struct
			x = sym('x',[1 2]);
			cofs.k = matlabFunction(symfun(self.k,x));
			cofs.r = matlabFunction(symfun(self.r,x));
			cofs.p = matlabFunction(symfun(self.p,x));

		end

	end

end
