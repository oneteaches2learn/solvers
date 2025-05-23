classdef ManufacturedFunctions2d_poisson < ManufacturedFunctions2d_elliptic

	properties
		r
		dr_du
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
			self@ManufacturedFunctions2d_elliptic(k,uTrue,NameValueArgs);

			% store additional coefficients
			self.r = r;

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

			% note: this is a kludge. The linear and nonlinear solvers make
			% different assumptions about how the reaction term is handled. So
			% that leads to different methods of manufacting the RHS. 
			if Coefficients.isNonlinear(self.r)
				f = self.divq + compose(self.r,self.uTrue);
			else
				f = self.divq + self.r * self.uTrue;
			end

		end
	end
end
