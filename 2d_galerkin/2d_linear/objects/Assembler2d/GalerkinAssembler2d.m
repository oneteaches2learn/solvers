classdef GalerkinAssembler2d

	properties
	end

	methods
		function self = GalerkinAssembler2d
		end
	end

	methods (Static)
		function dom = assembleDomainGeometry(xLim,yLim,inc,eps)

			% if two arguments, create unpunctured domain
			if nargin == 2
				dom = Domain2d(xLim,yLim);
			
			% if four arguments, create punctured domain
			elseif nargin == 4

				% if inclusion is passed as object
				if isa(inc,'Inclusion2d')

					% create punctured domain
					dom = Domain2d_punctured(xLim,yLim,inc,eps);

				% else if parameters for inclusion are passed as cell
				elseif isa(inc,'cell')

					% unpack parameters
					xLim_Y = inc{1};
					yLim_Y = inc{2};
					incRatio = inc{3};
					type = inc{4};

					% create inclusion object
					if strcmp(type,'circle')
						inc = Inclusion2d_circle(xLim_Y,yLim_Y,incRatio);
					elseif strcmp(type,'square')
						inc = Inclusion2d_square(xLim_Y,yLim_Y,incRatio);
					end

					% create punctured domain
					dom = Domain2d_punctured(xLim,yLim,inc,eps);
				end
			end
		end

		function dom = assembleBoundary(dom,bcTypes,varargin); 
		% assembleBoundary provides an interface with the Boundary2d object
		%  	Two input sets are allowed:
		%		(1) dom, bcTypes, bcConds
		%		(2) dom, bcTypes, u_D, u_N, alpha, u_R
		% In either case, bcTypes is a cell array, a string array, or a char
		% array. Each entry of bcTypes corresponds to an edge of the boundary.
		% The entries of bcTypes are 'D', 'N', 'R', or 'P', corresponding to
		% Dirichlet, Neumann, Robin, or Periodic boundary conditions,
		% respectively. For now, all boundaries will consist of four outer edges
		% and (possibly) some or many internal edges. If the same boundary
		% condition is to be applied to all internal edges, then the user can
		% pass an array of length 5, the first four entries corresponding to the
		% four internal edges, and the fifth entry being the boundary type
		% applied to all internal edges. 
		% 
		% The entries of bcConds are the boundary conditions themselves. In case
		% (1), the entries of bcConds are also listed by edge. Note that for
		% Robin boundary conditions, the corresponding input of bcConds should
		% be a nested cell array, with the first entry corresponding to the
		% coefficient function alpha and the second entry corresponding to u_R.
		% In case (2), the boundary conditions are given by type; in this case,
		% assembleBoundary sorts the boundary conditions by edge. In either
		% case, the inputs should be given as cell arrays. Entries of these cell
		% arrays may be of type double, sym, symfun, function_handle, or a
		% mixture of these. In any case, assembleBoundary will convert these
		% entries to function_handles. 
		
			% process bcTypes
			bcTypes = GalerkinAssembler2d.check_bcTypes(dom,bcTypes);
			
			% if inputs passed by boundary condition type, collate into bcConds
			if nargin == 6
				u_D 	= varargin{1};
				u_N 	= varargin{2};
				alpha 	= varargin{3};
				u_R 	= varargin{4};
				bcConds = GalerkinAssembler2d.get_bcConds(dom,bcTypes,u_D,u_N,alpha,u_R);
			
			% if inputs passed by BCtype, including dynamic BCs
			elseif nargin == 10
				u_D 	= varargin{1};
				u_N 	= varargin{2};
				alpha_R = varargin{3};
				u_R 	= varargin{4};
				alpha_T = varargin{5};
				beta_T 	= varargin{6};
				gamma_T = varargin{7};
				u_T 	= varargin{8};

				bcConds = GalerkinAssembler2d.get_bcConds( ...
							dom,bcTypes,u_D,u_N,alpha_R,u_R,alpha_T,beta_T,gamma_T,u_T);
			
			% if inputs passed by edge, store
			elseif nargin == 3
				bcConds = varargin{1};
			end

			% set boundary conditions on domain
			dom = dom.setBCTypes(bcTypes);
			dom = dom.setBCConditions(bcConds);

		end

		function dom = assembleNonlinearBoundary(dom,bcTypes,BCs)
		% assembleNonlinearBoundary adds d/du(BC) to the boundary conditions
		%
		% INPUTS
		%	dom: Domain2d object
		% 	bcTypes: char array with entries 'D', 'N', 'R', 'T', one for each
		% 		boundary edge
		%   BCs: cell array with symfun entries representing the boundary
		%   	conditions, one for each boundary edge
		%
		% The derivative of each boundary condition will be taken using symbolic
		% differentiation, so the inputs must be symbolic functions. For
		% dirichlet conditions, pass a 0 as the corresponding symbolic function.
		% Robin conditions are of the form:
		%
		%	du / dn = -alpha_R * (u - u_R).
		%
		% Entries of BCs corresponding to R condtion should be passed as a
		% nested cell array in the form {alpha_R,u_R}.
		
			% create necessary symbolic variables
			x = sym('x',[1 2],'real'); syms u; syms t;

			% converge doubles to symbolic functions
			for i = 1:length(BCs)

				% if R condition, extract u_R
				if bcTypes(i) == 'R'
					BCs{i} = BCs{i}{2};
				end
				
				% check if BC is a double
				if isnumeric(BCs{i})
					BCs{i} = symfun(BCs{i},x);
				end

				% take derivative of BC
				BCs{i} = symfun(diff(BCs{i},u),x);

				% assign to edge
				dom.boundary.edges(i).boundaryCondition_ddu = matlabFunction(BCs{i});

			end
		end
		
		function bcTypes = check_bcTypes(dom,bcTypes)
		% Let numEdges be the number of edges in the boundary of the domain.
		% Output bcTypes is a char array with length(bcTypes) = numEdges, each
		% char representing the BC type on the corresponding edge.
		% 
		% input bcTypes can be a character array or cell array of length 4, 5,
		% or length(bcTypes) = numEdges. Entries 1-4 correspond the four outer
		% edges of the domain. If the same BC is to be applied to all internal
		% edges, then the user can pass a length 5 array, the 5th entry
		% representing the BC type to apply to the internal edges. Otherwise,
		% length(bcTypes) = numEdges where each entry explicitly represents the
		% BC type on the corresponding edge. 
		%
		% check_bcTypes processes input bcTypes into output bcTypes per above.

			if iscell(bcTypes)
				for i = 1:length(bcTypes)
					temp(i) = bcTypes{i};
				end
				bcTypes = temp;
			end

			% extend bcTypes if necessary
			if length(bcTypes) ~= length(dom.boundary.edges)
				temp = repmat(bcTypes(5),1,length(dom.boundary.edges)-5);
				bcTypes = [bcTypes temp];
			end

		end

		function bcConds = get_bcConds(dom,bcTypes,u_D,u_N,alpha_R,u_R,alpha_T,beta_T,gamma_T,u_T)
		% get_bcConds sorts boundary conditions by edge
		%
		% Let numEdges be the number of edges in the boundary of the domain.
		%
		% bcTypes is a char array with length(bcTypes) = numEdges, each char
		% representing the BC type on the corresponding edge.
		%
		% each of u_D, u_N, alpha, u_R, beta, and u_T is a function_handle or a
		% cell array of length numEdges. For example, if the same function u_N =
		% g is to be applied to all Neumann-type edges in the domain, the user
		% can pass a single function handle. If different functions are to be
		% applied to different edges, the user can pass a cell array of length
		% numEdges, where each cell contains a function_handle representing the
		% boundary condition to be applied to the corresponding edge. 
		
			% collect boundary conditions as cell array
			for i = 1:length(bcTypes)

				% if edge i has Dirichlet BC
				if bcTypes(i) == 'D'
					if iscell(u_D), bcConds{i} = u_D{i}; 
					else, bcConds{i} = u_D; end

				% if edge i has Neumann BC
				elseif bcTypes(i) == 'N'
					if iscell(u_N), bcConds{i} = u_N{i};
					else, bcConds{i} = u_N; end

				% if edge i has Robin BC
				elseif bcTypes(i) == 'R'
					if iscell(alpha_R) && ~iscell(u_R)
						bcConds{i} = {alpha_R{i},u_R};
					elseif ~iscell(alpha_R) && iscell(u_R)
						bcConds{i} = {alpha_R,u_R{i}};
					elseif iscell(alpha_R) && iscell(u_R)
						bcConds{i} = {alpha_R{i},u_R{i}};
					else
						bcConds{i} = {alpha_R,u_R};
					end

				% if edge i has dynamic BC
				elseif bcTypes(i) == 'T'
					
					% fetch entries from cell arrays, if necessary
					if iscell(alpha_T), alpha_i = alpha_T{i};
					else alpha_i = alpha_T; end
					if iscell(beta_T), beta_i = beta_T{i};
					else beta_i = beta_T; end
					if iscell(gamma_T), gamma_i = gamma_T{i};
					else gamma_i = gamma_T; end
					if iscell(u_T), u_T_i = u_T{i};
					else u_T_i = u_T; end
					
					% store results
					bcConds{i} = {alpha_i,beta_i,gamma_i,u_T_i};
				end

			end
		end

		%{
		%NOTE: ORIGINAL VERSION
		function f = getFunctionHandles(f,vars)

			% convert to function handles
			if iscell(f)
				for i = 1:length(f)
					if iscell(f{i})
						func = f{i};
						alpha = func{1};
						u_R = func{2};
						f{i} = {matlabFunction(symfun(alpha,vars)),matlabFunction(symfun(u_R,vars))};
					else
						f{i} = matlabFunction(symfun(f{i},vars));
					end
				end
			elseif isnan(f)
				...
			else
				f = matlabFunction(symfun(f,vars));
			end

		end
		%}

		function f = getFunctionHandles(f,vars)
		% getFunctionHandles converts inputs f to function handles of the
		% variables vars. 
		%
		% INPUTS
		%	f: a function_handle, symfun, double, or cell array of these
		%	vars: a vector of symbolic variables, e.g. [x t] where x is a 2D
		%		vector and t is time.
		%
		% OUTPUTS
		%	f: a function_handle or cell array of function_handles
		%

			% convert to function handles
			if iscell(f)
				for i = 1:length(f)
					if iscell(f{i})
						func = f{i};
						alpha = func{1};
						u_R = func{2};
						f{i} = {GalerkinAssembler2d.convertToFunctionHandle(alpha,vars),...
								GalerkinAssembler2d.convertToFunctionHandle(u_R,vars)};
					else
						f{i} = GalerkinAssembler2d.convertToFunctionHandle(f{i},vars);
					end
				end
			elseif ~isa(f,'function_handle') && isnan(f)
				...
			else
				f = GalerkinAssembler2d.convertToFunctionHandle(f,vars);
			end

		end

		function f = convertToFunctionHandle(f,vars)

			if ~isa(f,'function_handle')
				f = matlabFunction(symfun(f,vars));
			end
		end

		function dom = assembleMesh(dom,p,base,NameValuePairs)

			arguments
				dom
				p
				base
				NameValuePairs.effectiveRegion = 'Omega_eps';
				NameValuePairs.meshInclusions = 'off';
			end

			% WHEN THE PROGRAM RANDOMLY FAILS, THIS IS THE PLACE WHERE THERE'S
			% PROBABLY AN ISSUE. CAN YOU PLEASE FIX THIS PERMANENTLY?
			if isa(dom,'Domain2d_punctured')
				dom = dom.setMesh(p,base, ...
					effectiveRegion = NameValuePairs.effectiveRegion, ...
					meshInclusions = NameValuePairs.meshInclusions);
			else
				dom = dom.setMesh(p,base);
			end

		end

		function auxfun = assembleCoefficients(k,r,f)
		% NOTE: you actually *want* to pass the coefficients as symfun functions
		% of x. Both elliptic and parabolic solvers are expecting
		% function_handles; constant coefficients will not work. Therefore,
		% constant coefficients must be converted to functions of space at
		% least. For parabolic problems in particular, the solver will check
		% which tensors need to be recomputed at each timestep, and this is done
		% by checking whether or not the coefficients used to construct the
		% tensor are / aren't time varying. If the coefficients are time varying
		% already, then they will remain time varying even after passing through
		% the code symfun(f,x). On the other hand, if the coefficient is
		% spatially varying only (or constant), then you do not want to impose
		% that it by time-varying, since the solver will then unnecessarily
		% recompute the corresponding tensor at each time step. Finally, if the
		% problem is elliptic, then there is no need for the coefficients to be
		% time varying. 
		
			% check function variables
			x = sym('x',[1 2],'real'); u = sym('u','real');

			% store diffusivity term
			if ~isa(k,'function_handle')
				k = symfun(k,x);
				auxfun.cofs.k = matlabFunction(k);
			else
				auxfun.cofs.k = k;
			end	
	
			% store reaction term
			if ~isa(r,'function_handle')
				r = symfun(r,x);
				auxfun.cofs.r = matlabFunction(r);
				auxfun.cofs.dr_du = matlabFunction(diff(r,u));
			else
				auxfun.cofs.r = r;
				auxfun.cofs.dr_du = 0;
				warning('Reaction term passed as function_handle. If your problem is to be solved using Newton iteration and r = r(u), then you must manually overwrite auxfun.cofs.dr_du.');
			end

			% store source term
			if ~isa(f,'function_handle')
				f = symfun(f,x);
				auxfun.f = matlabFunction(f);
			else
				auxfun.f = f;
			end

		end

		%{
		function f = assembleSource(f)

			% check function variables
			x = sym('x',[1 2],'real');
			f = symfun(f,x);

			% convert to function_handle
			f = matlabFunction(f);

		end
		%}
	end
end

