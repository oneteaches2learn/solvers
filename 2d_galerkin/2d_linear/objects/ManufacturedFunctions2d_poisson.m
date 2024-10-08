classdef ManufacturedFunctions2d_poisson < ManufacturedFunctions2d_elliptic

	properties
		r
		dr_du
	end

	methods
		function self = ManufacturedFunctions2d_poisson(k,r,uTrue)
			
			% call superclass constructor
			self@ManufacturedFunctions2d_elliptic(k,uTrue)

			% store additional coefficients
			self.r = r;

			if Coefficients.isNonlinear(r)
				u = sym('u','real');
				self.dr_du = diff(r,u);
			end

			% manufacture RHS
			self.f = self.manufactureRHS; 

		end

		function funcs = functionHandles(self)

			% call superclass method
			funcs = functionHandles@ManufacturedFunctions2d(self);

			% store additional functions
			x = sym('x',[1 2]);
			funcs.cofs.r = matlabFunction(symfun(self.r,x));

			% store derivatives of nonlinear functions
			if Coefficients.isNonlinear(self.r)
				funcs.cofs.dr_du = matlabFunction(symfun(self.dr_du,x));
			end

		end

		function f = manufactureRHS(self)
		% Manufactures symfun f (i.e. the right-hand side) from data

			% manufacture RHS
			x = sym('x',[1 2]); syms t;
			f = self.divq + self.r;

		end

	end

end
