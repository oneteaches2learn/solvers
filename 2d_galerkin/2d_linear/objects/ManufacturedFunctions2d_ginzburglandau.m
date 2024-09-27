classdef ManufacturedFunctions2d_ginzburglandau < ManufacturedFunctions2d_elliptic

	properties
		r
	end

	methods
		function self = ManufacturedFunctions2d_ginzburglandau(k,r,uTrue)
			
			% call superclass constructor
			self@ManufacturedFunctions2d_elliptic(k,uTrue)

			% store additional coefficients
			self.r = r;

			% manufacture RHS
			self.f = self.manufactureRHS; 

		end

		function funcs = functionHandles(self)

			% call superclass method
			funcs = functionHandles@ManufacturedFunctions2d(self);

			% store additional functions
			x = sym('x',[1 2]);
			funcs.cofs.r = matlabFunction(symfun(self.r,x));

		end

		function f = manufactureRHS(self)
		% Manufactures symfun f (i.e. the right-hand side) from data

			% manufacture RHS
			x = sym('x',[1 2]); syms t;
			f = self.divq - self.uTrue + self.uTrue^3;

		end

	end

end
