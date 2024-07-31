classdef ManufacturedFunctions2d_pennes < ManufacturedFunctions2d
% ManufacturedFunctions2d_pennes(p,k,r,uStar,uTrue) manufactures the source for a Pennes equation MMS test.
%
% The Pennes equation is 
%
%			(pu)_t - div (k grad u) + r(u - uStar) = f
%
% The ManufacturedFunctions2d_pennes subclass manufactures the source term f
% for this problem. 
%
% Author: Tyler Fara			Date: July 28, 2024
%-----------------------------------------------------------------------------%

	properties
		r	  % double or sym, coefficient in (pu)_t - div (k grad u) + r(u - uStar) = f
		uStar % double or sym, coefficient in (pu)_t - div (k grad u) + r(u - uStar) = f
	end

	methods
		function self = ManufacturedFunctions2d_pennes(p,k,r,uStar,uTrue)
		% ManufacturedFunctions2d_pennes(p,k,r,uStar,uTrue) inputs are symfun objects
			
			% call superclass constructor
			self@ManufacturedFunctions2d(p,k,uTrue)

			% store additional coefficients
			self.r = r;
			self.uStar = uStar;

			% manufacture RHS
			self.f = self.manufactureRHS; 

		end

		function f = manufactureRHS(self)
		% Manufactures symfun f (i.e. the right-hand side) from data

			% manufacture RHS
			x = sym('x',[1 2]); syms t;
			f = self.pu_t + self.divq + self.r * (self.uTrue - self.uStar);

		end

		function cofs = coefficients2FunctionHandles(self)
		% Converts symfun coefficients p, k, r, uStar to structure containing function_handles
		%	For cofs = self.coefficients2FunctionHandles outputs will be of
		%	format cofs.k, cofs.r, cofs.p, cofs.uStar;

			% convert coefficients to function handles, store as struct
			x = sym('x',[1 2]);
			cofs.k = matlabFunction(symfun(self.k,x));
			cofs.r = matlabFunction(symfun(self.r,x));
			cofs.p = matlabFunction(symfun(self.p,x));
			cofs.uStar = matlabFunction(symfun(self.uStar,x));

		end

	end

end
