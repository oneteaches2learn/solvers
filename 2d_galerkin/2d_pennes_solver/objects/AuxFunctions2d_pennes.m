classdef AuxFunctions2d_pennes < AuxFunctions2d

	properties
		uStar
	end

	methods
		function self = AuxFunctions2d_pennes(p,k,r,uStar,uInit,uTrue)
			
			% call superclass constructor
			self@AuxFunctions2d(p,k,r,uInit,uTrue)

			% store additional coefficients
			x = sym('x',[1 2]); syms t;
			self.uStar = symfun(uStar,[x t]);

			% manufacture RHS
			self.f = self.manufactureRHS; 

		end

		function f = manufactureRHS(self)

			% manufacture RHS
			x = sym('x',[1 2]); syms t;
			f = self.p * self.u_t + self.divq + self.r * (self.uTrue - self.uStar);

		end

		function cofs = coefficients2FunctionHandles(self)

			% convert coefficients to function handles, store as struct
			cofs.k = matlabFunction(self.k);
			cofs.r = matlabFunction(self.r);
			cofs.p = matlabFunction(self.p);
			cofs.uStar = matlabFunction(self.uStar);

		end

	end

end
