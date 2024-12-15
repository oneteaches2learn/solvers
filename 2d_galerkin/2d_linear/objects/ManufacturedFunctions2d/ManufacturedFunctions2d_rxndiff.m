classdef ManufacturedFunctions2d_rxndiff < ManufacturedFunctions2d_parabolic
% ManufacturedFunctions2d_rxndiff(p,k,r,uTrue) manufactures the source for a reaction-diffusion equation MMS test.
%
% The linear reaction-diffusion equation is 
%
%			(pu)_t - div (k grad u) + ru = f
%
% The ManufacturedFunctions2d_rxndiff subclass manufactures the source term f
% for this problem. 
%
% Author: Tyler Fara			Date: July 31, 2024
%-----------------------------------------------------------------------------%

	properties
		r	  % double or sym, coefficient in (pu)_t - div (k grad u) + ru = f
		dr_du
	end

	methods
		function self = ManufacturedFunctions2d_rxndiff(p,k,r,uTrue,NameValueArgs)
		% ManufacturedFunctions2d_pennes(p,k,r,uTrue) inputs are symfun objects
			
			arguments
				p
				k
				r
				uTrue
				NameValueArgs.u_N = 0;
				NameValueArgs.alpha_R = 0;
				NameValueArgs.u_R = 0;
			end
			
			% call superclass constructor
			self@ManufacturedFunctions2d_parabolic(p,k,uTrue,NameValueArgs);

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

		function f = manufactureRHS(self)
		% Manufactures symfun f (i.e. the right-hand side) from data

			% manufacture RHS
			x = sym('x',[1 2]); syms t;

			% note: this is a kludge. The linear and nonlinear solvers make
			% different assumptions about how the reaction term is handled. So
			% that leads to different methods of manufacting the RHS. 
			if Coefficients.isNonlinear(self.r)
				f = self.cu_t + self.divq + compose(self.r,self.uTrue);
			else
				f = self.cu_t + self.divq + self.r * self.uTrue;
			end


		end

		function funcs = functionHandles(self)

			% call superclass method
			funcs = functionHandles@ManufacturedFunctions2d_parabolic(self);

			% convert coefficients to function handles, store as struct
			x = sym('x',[1 2]);
			funcs.cofs.r = matlabFunction(symfun(self.r,x));

			% store derivatives of nonlinear functions
			funcs.cofs.dr_du = matlabFunction(symfun(self.dr_du,x));

		end

	end

end
