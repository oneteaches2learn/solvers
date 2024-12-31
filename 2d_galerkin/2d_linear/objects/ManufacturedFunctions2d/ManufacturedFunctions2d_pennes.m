classdef ManufacturedFunctions2d_pennes < ManufacturedFunctions2d_parabolic
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
		function self = ManufacturedFunctions2d_pennes(p,k,r,uStar,uTrue,NameValueArgs)
		% ManufacturedFunctions2d_pennes(p,k,r,uStar,uTrue) inputs are symfun objects
			
			arguments
				p
				k
				r
				uStar
				uTrue
				NameValueArgs.u_N = 0;
				NameValueArgs.alpha_R = 0;
				NameValueArgs.u_R = 0;
			end

			% call superclass constructor
			self@ManufacturedFunctions2d_parabolic(p,k,uTrue,NameValueArgs);

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
			f = self.cu_t + self.divq + self.r * (self.uTrue - self.uStar);

		end

		function funcs = functionHandles(self)

			% call superclass method
			funcs = functionHandles@ManufacturedFunctions2d_parabolic(self);

			% convert coefficients to function handles, store as struct
			x = sym('x',[1 2]);
			funcs.cofs.r = matlabFunction(symfun(self.r,x));
			funcs.cofs.uStar = matlabFunction(symfun(self.uStar,x));

		end

	end

end
