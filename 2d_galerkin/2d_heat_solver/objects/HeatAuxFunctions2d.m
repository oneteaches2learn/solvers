classdef HeatAuxFunctions2d

	properties
		p
		k
		r
		uInit
		q
		divq
		f
		uTrue
		coefficients
	end

	methods
		function self = HeatAuxFunctions2d(p,k,r,uInit,uTrue)
			
			% store inputs
			x = sym('x',[1 2]); syms t;
			self.p = symfun(p,[x t]);
			self.k = symfun(k,[x t]);
			self.r = symfun(r,[x t]);
			self.uTrue = symfun(uTrue,[x t]);
			self.uInit = symfun(self.uTrue(x(1),x(2),0),x);

			% manufacture flux
			q = -k * gradient(self.uTrue);
			q = formula(q);
			q = q(1:2);
			self.q = q;

			% manufacture divergence of flux
			q_x1 = diff(q(1),x(1));
			q_x2 = diff(q(2),x(2));
			divq = q_x1 + q_x2;
			self.divq = divq;

			% manufacture time derivative
			u_t = diff(uTrue,t);

			% manufacture RHS
			self.f = self.p * u_t + self.divq + self.r * self.uTrue;

		end

		function cofs = coefficients2FunctionHandles(self)

			% convert coefficients to function handles, store as struct
			cofs.k = matlabFunction(self.k);
			cofs.r = matlabFunction(self.r);
			cofs.p = matlabFunction(self.p);

		end

		function f = source2FunctionHandle(self)

			% convert source to function handle
			f = matlabFunction(self.f);

		end

		function uInit = uInit2FunctionHandle(self)

			% convert source to function handle
			uInit = matlabFunction(self.uInit);

		end

	end

end
