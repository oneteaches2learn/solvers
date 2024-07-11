classdef PoissonAuxFunctions2d

	properties
		k
		r
		q
		divq
		f
		uTrue
		coefficients
	end

	methods
		function self = PoissonAuxFunctions2d(k,r,uTrue)
			
			% store inputs
			x = sym('x',[1 2]);
			self.k = symfun(k,x);
			self.r = symfun(r,x);
			self.uTrue = symfun(uTrue,x);

			% manufacture flux
			q = -k * gradient(uTrue);
			q = formula(q);
			q = q(1:2);
			self.q = q;

			% manufacture divergence of flux
			q_x1 = diff(q(1),x(1));
			q_x2 = diff(q(2),x(2));
			divq = q_x1 + q_x2;
			self.divq = divq;

			% manufacture RHS
			self.f = self.divq + self.r * self.uTrue;

		end

		function cofs = coefficients2FunctionHandles(self)

			% convert coefficients to function handles, store as struct
			cofs.k = matlabFunction(self.k);
			cofs.r = matlabFunction(self.r);

		end

		function f = source2FunctionHandle(self)

			% convert source to function handle
			f = matlabFunction(self.f);

		end

	end

end
