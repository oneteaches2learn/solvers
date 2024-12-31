classdef ManufacturedFunctions2d < Coefficients
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
		k		% double or symfun, coefficient in (pu)_t - div (k grad u) = f
		q		% symfun, q = -k grad uTrue
		divq	% symfun, divq = div q = - div (k grad uTrue)
		f		% symfun, manufactured source function
		uTrue	% symfun, represents desired true solution
		u_N
		alpha_R
		u_R
		ODE
		bcConds
	end

	methods
		% CONSTRUCTOR
		function self = ManufacturedFunctions2d(k,uTrue,NameValueArgs)

			% store inputs
			x = sym('x',[1 2]);
			self.uTrue = symfun(uTrue,x);
			self.k = symfun(k,x);
			self.u_N = NameValueArgs.u_N;
			self.alpha_R = NameValueArgs.alpha_R;
			self.u_R = NameValueArgs.u_R;

			% manufacture data
			self.q = self.manufactureFlux;
			self.divq = self.manufactureFluxDivergence;

		end

		function q = manufactureFlux(self)
		% Manufactures symfun flux q = -k grad u

			% manufacture flux
			x = sym('x',[1 2]);
			uTrue = symfun(self.uTrue,x);
			q = -self.k * gradient(uTrue);
			q = formula(q);
			q = q(1:2);

		end

		function divq = manufactureFluxDivergence(self)
		% Manufactures symfun div q = - div (k grad u)

			% manufacture divergence of flux
			x = sym('x',[1 2]);
			q_x1 = diff(self.q(1),x(1));
			q_x2 = diff(self.q(2),x(2));
			divq = q_x1 + q_x2;

		end

		function f = manufactureRHS(self)
		% Manufactures symfun f (i.e. the right-hand side) from data

			% NOTE: placeholder function. Actually handled by specific subclasses.
			...

		end

		function funcs = functionHandles(self)

			% add outputs
			x = sym('x',[1 2]);
			funcs.cofs.k = matlabFunction(symfun(self.k,x));
			funcs.f = matlabFunction(symfun(self.f,x));

		end

		function self = manufactureBCs(self,dom,bcTypes,auxfun,varargin);

			% unpack coefficients
			uTrue = self.uTrue;
			%if isa(self,'GalerkinMMS2d_parabolic')
			if 1 == 1
				u_t = self.u_t;
			end
			q = self.q;
			nEdges = dom.boundary.nEdges;

			% set symbolic variables
			%if isa(self,'GalerkinMMS2d_parabolic')
			if 1 == 1
				x = sym('x',[1 2]); syms t;
				vars = [x t];
			else
				x = sym('x',[1 2]);
				vars = x;
			end

			% set edge normal vectors
			dom = self.setEdgeNormalVectors_outerBoundary(dom);
			if isa(dom,'Domain2d_punctured')
				dom = self.setEdgeNormalVectors_inclusions(dom);
			end

			% extend bcTypes
			bcTypes_interior = bcTypes(1:4);
			bcTypes_exterior = bcTypes(5:end);
			if length(bcTypes_exterior) == 1
				bcTypes_exterior = repmat(bcTypes_exterior,1,nEdges-4);
			end
			bcTypes = [bcTypes_interior, bcTypes_exterior];

			% Manufacture Dirichlet BC
			u_d = uTrue;
			%u_d = matlabFunction(u_d);

			% assign boundary functions
			for i = 1:nEdges

				% assign Dirichlet BC
				if strcmp(bcTypes(i),'D')
					bcConds{i} = u_d;

				% assign Neumann BC
				elseif strcmp(bcTypes(i),'N')
					n_i = symfun(dom.boundary.edges(i).outwardNormal,vars);
					g_i = symfun(sum(q.*n_i),vars);
					bcConds{i} = g_i;

				% assign Robin BC
				elseif strcmp(bcTypes(i),'R')
					alpha_i = symfun(1.0,vars);
					n_i = symfun(dom.boundary.edges(i).outwardNormal,vars);
					g_i = symfun(uTrue - sum(q .* n_i) / alpha_i,vars);
					bcConds{i} = {alpha_i,g_i};

				%{
				% assign dynamic BC
				elseif dom.boundary.edges(i).boundaryType == 'T'
					alpha_i = symfun(1.0,vars);
					beta_i  = symfun(1.0,vars);
					gamma_i = symfun(1.0,vars);
					n_i = symfun(dom.boundary.edges(i).outwardNormal,vars);
					g_i = symfun((alpha_i * uTrue + beta_i * u_t - sum(q .* n_i)) / gamma_i,vars);
					alpha_i = matlabFunction(alpha_i);
					beta_i  = matlabFunction(beta_i);
					gamma_i = matlabFunction(gamma_i);
					g_i = matlabFunction(g_i);
					dom.boundary.edges(i).boundaryCondition = {alpha_i,beta_i,gamma_i,g_i};
				%}
				end
			end

			self.bcConds = bcConds;

		end

		function dom = setEdgeNormalVectors_outerBoundary(self,dom)

			% store variables
			dl = dom.boundary.dl.mat;

			% get normal vectors
			n_lower = [0; -1];
			n_right = [1; 0];
			n_upper = [0; 1];
			n_left  = [-1; 0];

			n_vectors = [n_lower, n_right, n_upper, n_left];

			% loop over columns of decomposed geometry description matrix
			for j = 1:4;

				% for now, it is assumed the outer boundary is rectangular
				if dl(1,j) == 2

					n = n_vectors(:,mod(j-1,4)+1);

				end
				
				% store normal vector
				dom.boundary.edges(j).outwardNormal = n;

			end
		end

		function dom = setEdgeNormalVectors_inclusions(self,dom)

			% store variables
			dl = dom.boundary.dl;
			edgeIDs = dl.segEdgeDict;
			scale_eps = 1 / dom.epsilon;

			% setup symbolic functions for circle edges
			x = sym('x',[1 2],'real');

			% if unit vectors to circle are needed, store in advance
			if sum(find(dl.mat(1,:) == 1)) > 0
				n_lower = dom.inclusion.Q.unitNormal_lower;
				n_upper = dom.inclusion.Q.unitNormal_upper;
			end

			for i = 1:length(dl.edgeSegDict_inclusions)

				% get first segment of edge
				seg = dl.edgeSegDict_inclusions{i}(1);

				%if segment is of type: circle
				if dl.segType(seg) == 1

					% set x-translation
					temp = dl.circleCenter(seg);
					x_translate = temp(1);
					x_transformed = scale_eps * (x(1) - x_translate);

					% if lower edge of circle
					if dl.isLowerCircle(i)

						n = symfun(n_lower(x_transformed,x(2)),x);

					% else upper edge of circle
					else

						n = symfun(n_upper(x_transformed,x(2)),x);

					end

				% if segment is of type: line
				elseif dl.segType(seg) == 2

					% for now, presume there are four outer edges
					edgeNum = i + 4;
					orientation = dl.squareEdgeOrientation(edgeNum);
					n = dom.inclusion.Q.unitNormal(:,orientation);
					
				end
				
				% store normal vector
				% for now, presume there are four outer edges
				dom.boundary.edges(i + 4).outwardNormal = n;

			end

		end

		function self = set.ODE(self,ode)
			self.ODE = ode;
		end

	end

end
