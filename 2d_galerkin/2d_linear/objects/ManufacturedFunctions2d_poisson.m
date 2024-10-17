classdef ManufacturedFunctions2d_poisson < ManufacturedFunctions2d_elliptic

	properties
		r
		dr_du
		u_N
		alpha_R
		u_R
	end

	methods
		function self = ManufacturedFunctions2d_poisson(k,r,uTrue,NameValueArgs)

			arguments
				k
				r
				uTrue
				NameValueArgs.u_N = 0;
				NameValueArgs.alpha_R = 0;
				NameValueArgs.u_R = 0;
			end
			
			% call superclass constructor
			self@ManufacturedFunctions2d_elliptic(k,uTrue)

			% store additional coefficients
			self.r = r;
			self.u_N = NameValueArgs.u_N;
			self.alpha_R = NameValueArgs.alpha_R;
			self.u_R = NameValueArgs.u_R;

			% compute derivatives of nonlinear functions
			if Coefficients.isNonlinear(r)
				u = sym('u','real');
				self.dr_du = diff(r,u);
			else
				self.dr_du = 0;
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
			funcs.cofs.dr_du = matlabFunction(symfun(self.dr_du,x));

		end

		function f = manufactureRHS(self)
		% Manufactures symfun f (i.e. the right-hand side) from data

			% manufacture RHS
			x = sym('x',[1 2]); syms t;
			f = self.divq + compose(self.r,self.uTrue);

		end

	end

end
