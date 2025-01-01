classdef Coefficients

	methods
		function self = Coefficients
		end
	end

	methods (Static)
		function [f,t,u] = checkVariables(f)
		% checkVariables(f) converts f to an anonymous function handle in the
		% variables x, t, and u.
		% 
		% INPUTS
		%	f - double, sym or symfun, or function handle
		%
		% OUTPUTS
		%	f - function handle, f(x1,x2,t,u)
		%	t - double, time variable
		%	u - double, nonlinear variable
		%
		% In this package, the most general need is for f to be a nonlinear
		% function of space and time, where space is two dimensional. To support
		% this, all function handles are presumed to be functions of x, t, and
		% u. This assumption is made even when solving linear problems, or
		% stationary problems. In the case that f is constant, linear, or
		% stationary, the additional variables are added and set to zero. In the
		% case that f is not an anonymous function handle, i.e. if f is a
		% double, sym, or symfun, then f is first converted to a function
		% handle.

			if ~Coefficients.isNonlinear(f) && ~Coefficients.isTimeVarying(f)
				f = @(x1,x2,t,u)(f(x1,x2));
				t = 0;
				u = 0;
			elseif ~Coefficients.isNonlinear(f) && Coefficients.isTimeVarying(f)
				f = @(x1,x2,t,u)(f(x1,x2,t));
				u = 0;
			elseif Coefficients.isNonlinear(f) && ~Coefficients.isTimeVarying(f)
				f = @(x1,x2,t,u)(f(x1,x2,u));
				t = 0;
			end
		end

		function f = wrap2d(f)
		% wrap2d(f) ensures that f returns an output of the same size as the
		% input x1.
		%
		% INPUTS
		%	f - function handle, f(x1,x2,t,u)
		%
		% OUTPUTS
		%	f - function handle, f(x1,x2,t,u) + zeros(size(x1))
		%
		% In this package, f will usually be used to return a function at each
		% node of some mesh; i.e. the output of f will usually be a matrix.
		% However, if f is constant in the variable x, then MATLAB's
		% broadcasting rules will cause f to return a scalar. To ensure that f
		% returns a matrix, wrap2d adds an appropriately sized matrix of zeros
		% to the output of f.
		
			f = @(x1,x2,t,u)(f(x1,x2,t,u) + zeros(size(x1)));
		end


		function result = isTimeVarying(f)
		% isTimeVarying(f) returns true if the function handle f is time varying.
		%
		% INPUTS
		%	f - double, sym or symfun, or function handle
		%
		% OUTPUTS
		%	result - logical, true if f is time varying, false otherwise
		%
		% isTimeVarying checks if input f is a function of variable t. To do so,
		% isTimeVarying must be converted into an anonymous function handle. If f
		% is a double, then f represents a constant function and is not time
		% varying. In this case, f is converted to a symfun in the spatial
		% variable x. If f is a sym or symfun, then f is converted to a function
		% handle. The function handle is then checked for the presence of the
		% variable t.

			if isa(f,'double')
				x = sym('x',[1 2]);
				f = symfun(f,x);
			end
			
			if isa(f,'sym') || isa(f,'symfun')
				x = sym('x',[1 2]);
				f = matlabFunction(f,x);
			end

			result = Coefficients.hasVariable(f,'t');
		end

		function result = isNonlinear(f)
		% isNonlinear(f) returns true if the function handle f is nonlinear in
		% the variable u.
		%
		% INPUTS
		%	f - double, sym or symfun, or function handle
		%
		% OUTPUTS
		%	result - logical, true if f is nonlinear in u, false otherwise
		%
		% isNonlinear checks if input f is a function of variable u. To do so,
		% isNonlinear must be converted into an anonymous function handle. If f
		% is a double, then f represents a constant function and is nonlinear in
		% u. In this case, f is converted to a symfun in the spatial variable x.
		% If f is a sym or symfun, then f is converted to a function handle. The
		% function handle is then checked for the presence of the variable u.
		
			if isa(f,'double')
				x = sym('x',[1 2]);
				f = symfun(f,x);
			end

			if isa(f,'sym') || isa(f,'symfun')
				f = matlabFunction(f);
			end
			
			result = Coefficients.hasVariable(f,'u');
		end

		function result = isCoupled(f)
		% isNonlinear(f) returns true if the function handle f is a function of
		% the ODE variable v.
		%
		% INPUTS
		%	f - double, sym or symfun, or function handle
		%
		% OUTPUTS
		%	result - logical, true if f is a function of v, false otherwise
		%
		% isNonlinear checks if input f is a function of variable v. To do so,
		% isNonlinear must be converted into an anonymous function handle. If f
		% is a double, then f represents a constant function and is not a
		% function of v. In this case, f is converted to a symfun in the spatial
		% variable x.  If f is a sym or symfun, then f is converted to a
		% function handle. The function handle is then checked for the presence
		% of the variable v.
		
			if isa(f,'double')
				x = sym('x',[1 2]);
				f = symfun(f,x);
			end

			if isa(f,'sym') || isa(f,'symfun')
				f = matlabFunction(f);
			end
			
			result = Coefficients.hasVariable(f,'v');
		end

		function result = hasVariable(f,checkVar)
		% hasVariable(f,checkVar) returns true if the function handle f is a
		% function of the variable checkVar.
		%
		% INPUTS
		%	f - double, sym or symfun, or function handle
		%	checkVar - string, variable to check for
		%
		% OUTPUTS
		%	result - logical, true if f is a function of checkVar, false otherwise
		%
		% In this package, functions are presumed to be functions of x, t, and
		% u. checkVar could, presumably, be some other variable. But hasVariable
		% has been designed and tested for checkVar = 't' or 'u'.
		
			if isa(f,'double')
				x = sym('x',[1 2]);
				f = symfun(f,x);
			end
			
			if isa(f,'sym') || isa(f,'symfun')
				f = matlabFunction(f,x);
			end
			
			usedVars = Coefficients.getFunctionVariables(f);
			result = any(strcmp(usedVars,checkVar));
		end

		function vars = getFunctionVariables(fh)
		% getFunctionVariables(fh) returns the variables of function handle fh
		%
		% INPUTS
		%	fh - function handle
		%
		% OUTPUTS
		%	vars - cell array of strings, variables of fh
		%
		% getFunctionVariables uses regex to collect variable names from
		% function_handle fh. For example, if fh is a function of x, t, and u,
		% then vars will be {'x1', 'x2', 't','u'}. Presumably,
		% getFunctionVariables could check for variables other than x, t, and u.
		% But getFunctionVariables has been designed and tested for x, t, and u.
		
			% use regex to collect variable names from function_handle
			code = func2str(fh);
			pattern = '@\((.*?)\)';
			tokens = regexp(code,pattern,'tokens');

			if ~isempty(tokens)
				vars = strsplit(tokens{1}{1},',');
				vars = strtrim(vars);

			else
				vars = {};

			end
		end

	end

end

