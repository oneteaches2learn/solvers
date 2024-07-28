classdef AuxFunctions2d

	properties
		p
		k
		uInit
		q
		divq
		u_t
		f
		uTrue
		coefficients
	end

	methods
		function self = AuxFunctions2d(p,k,uInit,uTrue)
			
			% store inputs
			x = sym('x',[1 2]); syms t;
			self.p = symfun(p,[x t]);
			self.k = symfun(k,[x t]);
			self.uTrue = symfun(uTrue,[x t]);
			self.uInit = symfun(self.uTrue(x(1),x(2),0),x);

			% manufacture data
			self.q = self.manufactureFlux;
			self.divq = self.manufactureFluxDivergence;
			self.u_t = self.manufactureTimeDerivative;

		end

		function cofs = coefficients2FunctionHandles(self)

			% convert coefficients to function handles, store as struct
			cofs.k = matlabFunction(self.k);
			cofs.p = matlabFunction(self.p);

		end

		function q = manufactureFlux(self)

			% manufacture flux
			q = -self.k * gradient(self.uTrue);
			q = formula(q);
			q = q(1:2);

		end

		function divq = manufactureFluxDivergence(self)

			% manufacture divergence of flux
			x = sym('x',[1 2]); syms t;
			q_x1 = diff(self.q(1),x(1));
			q_x2 = diff(self.q(2),x(2));
			divq = q_x1 + q_x2;

		end

		function u_t = manufactureTimeDerivative(self)

			% manufacture time derivative
			x = sym('x',[1 2]); syms t;
			u_t = diff(self.uTrue,t);

		end

		function f = source2FunctionHandle(self)

			% convert source to function handle
			f = matlabFunction(self.f);

		end

		function uInit = uInit2FunctionHandle(self)

			% convert uInit to function handle
			uInit = matlabFunction(self.uInit);

		end

		function uTrue = uTrue2FunctionHandle(self)

			% convert uTrue to function handle
			uTrue = matlabFunction(self.uTrue);

		end

	end

end
