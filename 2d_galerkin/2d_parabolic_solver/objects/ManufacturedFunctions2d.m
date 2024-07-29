classdef ManufacturedFunctions2d
% ManufacturedFunctions2d(p,k,utrue) manufactures the source for a parabolic MMS test.
%
% The simplest parabolic problem is of the form
%
%			(pu)_t - div (k grad u) = f
%
% In an MMS test, uTrue and the coefficients k and p are presumed known.
% Source f is 'manufactured' from these quantities. The ManufacturedFunctions2d
% class manufactures and stores f, as well as flux q = - k div u, the
% divergence of flux divq, and the time derivative u_t.
%
% Author: Tyler Fara			Date: July 28, 2024
%-----------------------------------------------------------------------------%
% Notes
%  (1) ManufacturedFunctions2d takes advantage of MATLAB's symbolic functions
%  toolkit to compute derivatives and to assemble the resulting source function
%  f. As a result, the inputs to the ManufacturedFunctions2d constructor and
%  all of the properties of ManufacturedFunctions2d are symfun objects.
%
%  (2) Ultimately, the manufactured source f and the coefficients p and k will
%  be given to a solver, and solvers take function_handle objects. Thus, the
%  ManufacturedFunctions2d class contains methods to convert each symfun to its
%  corresponding function_handle. 
%
%  (3) ManufacturedFunctions2d is a superclass that is not directly used by any
%  solver. Instead, subclasses of ManufacturedFunctions2d exist that are
%  tailor-made for each PDE. 
%
%  (4) The f = self.manufactureRHS method is a placeholder method. Each
%  subclass of ManufacturedFunctions2d will have a specific manufactureRHS
%  method that is tailor-made for the PDE that the subclass is designed to
%  handle. 
%
%  (5) The coefficients2FunctionHandles method is another placeholder method.
%  Each subclass of ManufacturedFunctions2d will have a specific
%  coefficients2FunctionHandles method that is, again, tailor-made for the
%  problem at hand. 
%
%  (6) Also, the coefficients2FunctionHandles method is is one of many methods
%  that converts symfuns to function_handles. However, unlike all of the other
%  such methods, which all output a function_handle object, the
%  coefficients2FunctionHandles method outputs a structure, each property of
%  which is a function_handle. This was done because one cannot know a-priori
%  how many coefficients a parabolic PDE will have, though it will certainly
%  have at least two coefficients. (I.e., the simplest parabolic PDE, the heat
%  equation, listed above, already has two coefficients p and k.) Outputting a
%  structure allows flexibility for each subclass to load as many coefficients
%  as is necessary into the structure.  
%-----------------------------------------------------------------------------%


	properties
		p		% double or symfun, coefficient in (pu)_t - div (k grad u) = f
		k		% double or symfun, coefficient in (pu)_t - div (k grad u) = f
		uInit	% symfun, uInit(x) = uTrue(x,0)
		q		% symfun, q = -k grad uTrue
		divq	% symfun, divq = div q = - div (k grad uTrue)
		u_t		% symfun, u_t is time derivative of uTrue
		pu_t	% symfun, u_t is time derivative of p * uTrue
		f		% symfun, manufactured source function
		uTrue	% symfun, represents desired true solution
	end

	methods
		function self = ManufacturedFunctions2d(p,k,uTrue)
		% ManufacturedFunctions2d(p,k,uTrue) inputs are symfun objects
			
			% store inputs
			x = sym('x',[1 2]); syms t;
			self.p = p;
			self.k = k;
			self.uTrue = symfun(uTrue,[x t]);
			self.uInit = symfun(self.uTrue(x(1),x(2),0),x);

			% manufacture data
			self.q = self.manufactureFlux;
			self.divq = self.manufactureFluxDivergence;
			self.u_t = self.manufactureTimeDerivative;
			self.pu_t = self.manufactureTimeDerivative2;

			% manufacture RHS: handled by subclasses
				... 

		end


		function q = manufactureFlux(self)
		% Manufactures symfun flux q = -k grad u

			% manufacture flux
			q = -self.k * gradient(self.uTrue);
			q = formula(q);
			q = q(1:2);

		end

		function divq = manufactureFluxDivergence(self)
		% Manufactures symfun div q = - div (k grad u)

			% manufacture divergence of flux
			x = sym('x',[1 2]); syms t;
			q_x1 = diff(self.q(1),x(1));
			q_x2 = diff(self.q(2),x(2));
			divq = q_x1 + q_x2;

		end

		function u_t = manufactureTimeDerivative(self)
		% Manufactures symfun u_t, the time derivative of u

			% manufacture time derivative
			x = sym('x',[1 2]); syms t;
			u_t = diff(self.uTrue,t);

		end

		function pu_t = manufactureTimeDerivative2(self)
		% Manufactures symfun u_t, the time derivative of u

			% manufacture time derivative
			x = sym('x',[1 2]); syms t;
			pu_t = diff(self.p * self.uTrue,t);

		end

		function f = manufactureRHS(self)
		% Manufactures symfun f (i.e. the right-hand side) from data

			% NOTE: placeholder function. Actually handled by specific subclasses.
			...

		end

		function f = manufactureRHS2(self)
		% Manufactures symfun f (i.e. the right-hand side) from data

			% NOTE: placeholder function. Actually handled by specific subclasses.
			...

		end

		function cofs = coefficients2FunctionHandles(self)
		% Converts symfun coefficients p and k to structure containing function_handles

			% NOTE: placeholder function. Actually handled by specific subclasses.
			...

		end

		function f = source2FunctionHandle(self)
		% Converts symfun f to function_handle

			% convert source to function handle
			f = matlabFunction(self.f);

		end

		function uInit = uInit2FunctionHandle(self)
		% Converts symfun uInit = u(x,0) to function_handle

			% convert uInit to function handle
			uInit = matlabFunction(self.uInit);

		end

		function uTrue = uTrue2FunctionHandle(self)
		% Converts symfun uTrue to function_handle

			% convert uTrue to function handle
			uTrue = matlabFunction(self.uTrue);

		end

	end

end
