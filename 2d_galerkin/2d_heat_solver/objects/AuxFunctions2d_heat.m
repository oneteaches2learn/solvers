classdef AuxFunctions2d_heat < AuxFunctions2d

	properties
	end

	methods
		function self = AuxFunctions2d_heat(p,k,r,uInit,uTrue)
			
			% call superclass constructor
			self@AuxFunctions2d(p,k,r,uInit,uTrue)

			% manufacture RHS
			self.f = self.manufactureRHS; 

		end

		function f = manufactureRHS(self)

			% manufacture RHS
			x = sym('x',[1 2]); syms t;
			f = self.p * self.u_t + self.divq + self.r * self.uTrue;

		end


		function cofs = coefficients2FunctionHandles(self)

			% convert coefficients to function handles, store as struct
			cofs.k = matlabFunction(self.k);
			cofs.r = matlabFunction(self.r);
			cofs.p = matlabFunction(self.p);

		end

	end

end
