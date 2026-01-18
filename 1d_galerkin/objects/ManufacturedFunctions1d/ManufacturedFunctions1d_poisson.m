classdef ManufacturedFunctions1d_poisson < ManufacturedFunctions1d_elliptic

	properties
		r
		dr_du
		BCtypes
	end

	methods
		function self = ManufacturedFunctions1d_poisson(k,r,uTrue,BCtypes,NameValueArgs)

			arguments
				k
				r
				uTrue
				BCtypes = '';
				NameValueArgs.u_N = 0;
				NameValueArgs.alpha_R = 0;
				NameValueArgs.u_R = 0;
				NameValueArgs.reaction_form = 'mass'
			end
			
			% call superclass constructor
			self@ManufacturedFunctions1d_elliptic(k,uTrue);

			% store additional coefficients
			self.r = r;

			% compute derivatives of nonlinear functions
			if Coefficients.isNonlinear(r)
				u = sym('u','real');
				self.dr_du = diff(r,u);
			else
				self.dr_du = 0;
			end

			% if boundary conditions are provided, store them
			if ~isempty(BCtypes)
				self.BCtypes = BCtypes;
			end

			% store reaction form
			self.reaction_form = NameValueArgs.reaction_form;

			% manufacture RHS
			self.f = self.manufactureRHS; 

		end


		% MANUFACTURING FUNCTIONS
		function f = manufactureRHS(self)
		% Manufactures symfun f (i.e. the right-hand side) from data
		%
		% NOTE: ManufacturedFunctions1d manufactures f for either of
		% 
		% 	(1) - div (k grad u) + r(u) = f  (nonlinear, 'source' form)
		% 	(2) - div (k grad u) + r*u = f   (linear, 'mass' form)
		%
		% so this method must account for both cases. It does so by checking
		% whether r = r(u), in which case it assumes case (1). Or, if the
		% property reaction_form has been manually overwritten as reaction_form
		% = 'source', then ManufacturedFunctions1d_poisson will also assume case
		% (1). Otherwise, it assumes case (2). 
		
			% manufacture RHS
			syms x; syms t; syms u;

			if Coefficients.isNonlinear(self.r) || strcmp(self.reaction_form, "source")
				%f = self.divq + compose(self.r,self.uTrue);
				f = self.divq + compose(self.r,self.uTrue,u,x);
			else
				f = self.divq + self.r * self.uTrue;
			end

		end


		% OUTPUTTER FUNCTIONS
		function funcs = functionHandles(self)

			% call superclass method
			funcs = functionHandles@ManufacturedFunctions1d(self);

			% store additional functions
			syms x;
			funcs.cofs.r = matlabFunction(symfun(self.r,x));

			% store derivatives of nonlinear functions
			funcs.cofs.dr_du = matlabFunction(symfun(self.dr_du,x));

		end

		function funcs = outputCoefficients(self)
		% Alias for functionHandles method

			funcs = self.functionHandles;

		end

		function BCs = outputBoundaryConditions(self,dom,BCtypes)

			% store variables
			% input handling
			if nargin < 3
				BCtypes = self.BCtypes;
			end

			% loop over boundary types
			for i = 1:length(BCtypes)

				% assign Dirichlet BC
				if BCtypes(i) == 'D'
					BCvals{i} = matlabFunction(symfun(self.uTrue, sym('x','real')));
				end

				% assign Neumann BC
				if BCtypes(i) == 'N'
					% determine outward normal
					if i == 1, n = -1; else, n = 1; end

					% compute Neumann BC value
					bc = n * self.q;
					BCvals{i} = matlabFunction(symfun(bc, sym('x','real')));
				end
					
			end

			% create Boundary1d object
			BCs = Boundary1d(dom.xLim,BCtypes,BCvals);

			% for manufactured BCs, d/du is zero
			BCs.bConditions_du = {0,0};

		end

	end
end
