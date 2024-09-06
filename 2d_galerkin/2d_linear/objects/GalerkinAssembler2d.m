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
				u_D = varargin{1};
				u_N = varargin{2};
				alpha = varargin{3};
				u_R = varargin{4};
				bcConds = GalerkinAssembler2d.get_bcConds(dom,bcTypes,u_D,u_N,alpha,u_R);
			
			% if inputs passed by edge, store
			elseif nargin == 3
				bcConds = varargin{1};
			end

			% set boundary conditions on domain
			dom = dom.setBCTypes(bcTypes);
			dom = dom.setBCConditions(bcConds);

		end

		function bcTypes = check_bcTypes(dom,bcTypes)

			% convert cell to char array, if necessary
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

		function bcConds = get_bcConds(dom,bcTypes,u_D,u_N,alpha,u_R)
		% get_bcConds sorts boundary conditions by edge
		
			% collect boundary conditions as cell array
			for i = 1:length(bcTypes)
				if bcTypes(i) == 'D'
					if iscell(u_D), bcConds{i} = u_D{i}; 
					else, bcConds{i} = u_D; end
				elseif bcTypes(i) == 'N'
					if iscell(u_N), bcConds{i} = u_N{i};
					else, bcConds{i} = u_N; end
				elseif bcTypes(i) == 'R'
					if iscell(alpha) && ~iscell(u_R)
						bcConds{i} = {alpha{i},u_R};
					elseif ~iscell(alpha) && iscell(u_R)
						bcConds{i} = {alpha,u_R{i}};
					elseif iscell(alpha) && iscell(u_R)
						bcConds{i} = {alpha{i},u_R{i}};
					else
						bcConds{i} = {alpha,u_R};
					end
				end
			end
		end

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

		function dom = assembleMesh(dom,p,base)

			% set mesh on domain
			dom = dom.setMesh(p,base);

		end

		function cofs = assembleCoefficients(p,k)

			% check function variables
			x = sym('x',[1 2],'real'); syms t;
			p = symfun(p,x);
			k = symfun(k,x);

			% convert to function_handles
			cofs.p = matlabFunction(p);
			cofs.k = matlabFunction(k);

		end

		function uInit = assembleInitialCondition(uInit)

			% check function variables
			x = sym('x',[1 2],'real'); syms t;
			uInit = symfun(uInit,[x t]);

			% convert to function_handles
			uInit = matlabFunction(uInit);

		end

		function time = assembleTimeStepping(T,dt,eq)

			% create time-stepping object
			if nargin == 2
				time = TimeStepping(T,dt);
			elseif nargin == 3
				time = TimeStepping(T,dt,eq);
			end

			% set time-stepping mesh
			time = time.setMesh;

		end

		function f = assembleSource(f)

			% check function variables
			x = sym('x',[1 2],'real'); syms t;
			f = symfun(f,x);

			% convert to function_handle
			f = matlabFunction(f);

		end
	end
end
