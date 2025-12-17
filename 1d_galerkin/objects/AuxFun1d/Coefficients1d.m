classdef Coefficients1d < Coefficients

    properties
    end
        
    methods
        function self = Coefficients1d()
            % Constructor (empty)
        end
    end 

    methods (Static)        
        function func = wrapFunction(func)
        % WRAPFUNCTION Ensure that the input function is a function handle.
        % USAGE:
        %   func = wrapFunction(self,func)
        % INPUTS:
        %   func - function to be wrapped (can be function handle, symfun, sym, or constant)
        % OUTPUTS:
        %   func - wrapped function handle  
        %
        % NOTE:
        %   (1) Wrapping the function is necessary in the case where the input is
        %   constant either in time or in space. The expected output of any
        %   function evaluation is a matrix of size nSpaceEvals x nTimeEvals.
        %   But if a function is constant in space or time, MATLAB will return a
        %   single value an array of constants. Wrapping the function ensures that
        %   the output is always of the correct size. 
        %
        %   (2) A function is `wrapped' by adding an array of zeros of the same
        %   size as the expected output. However, the exact dimensions of this
        %   output array do not need to be known at the time of wrapping.
        %   Rather, by adding: 
        %       
        %       func = @(x,t) func(x,t) + zeros(x,t);
        %   
        %   a zero array of the correct size will be created at the time of
        %   function evaluation. This is sufficient to ensure that the output is
        %   always of the correct size.
        %
        %   (3) The expected input of any function evaluation is a
        %   function_handle, so wrapFunction also converts symfun, sym, or
        %   constant (i.e. double) inputs to a function handle. In the case of a
        %   constant function, evauating a function_handle is slightly more
        %   expensive than just creating an array of the appropriate size.
        %   However, this extra computational cost is worth it for the sake of
        %   reducing code complexity. 
        
            %{
            % if func is a constant (i.e. double) convert to function_handle
            if isa(func,'double') && length(func) == 1
                func = @(x,t) func;
            elseif isa(func,'double') && length(func) > 1
                error('If func is a double, it must be a scalar value.');
            end

            % if func is symfun convert to function_handle
            if isa(func,'symfun')
                syms x t
                func = matlabFunction(symfun(func,[x t]));
            end

            % if func is sym convert to function_handle
            if isa(func,'sym')
                syms x t
                func = matlabFunction(symfun(func,[x t]));
            end

            % if func is function_handle but only of x, convert to function_handle of (x,t)
            if isa(func,'function_handle') 
                func = Domain1d.checkVariables(func);
            end
            %}
            %{
            % convert to function handle if necessary
            func = Coefficients1d.func2functionHandle(func);

            % check variable dependence and update function signature
            func = Coefficients1d.checkVariables(func);
            %}

            % wrap function handle
            func = @(x,t) func(x,t) + zeros(length(x),length(t));

        end

        function func = func2functionHandle(func)
        % FUNC2FUNCTIONHANDLE Convert input function to function handle.
        % USAGE:
        %   func = func2functionHandle(func)
        % INPUTS:
        %   func - function to be converted (can be function handle, symfun, sym, or constant)
        % OUTPUTS:
        %   func - converted function handle

            % func is symfun 
            if isa(func,'symfun')
                syms x
                func = matlabFunction(symfun(func,x));

            % func is sym 
            elseif isa(func,'sym')
                syms x
                func = matlabFunction(symfun(func,x));

            % func is length 1 double with Symbolic Toolbox
            elseif isnumeric(func) ...
                && length(func) == 1  ...
                && license('test', 'Symbolic_Toolbox')

                syms x;
                func = matlabFunction(symfun(func,x));

            % func is length 1 double without Symbolic Toolbox
            elseif isnumeric(func) ...
                && length(func) ==1 ...
                && ~license('test', 'Symbolic_Toolbox')

                func = @(x) func;
            end

            % check that func is a function_handle
            if ~isa(func,'function_handle')
                error('Unable to convert input to function handle.');
            end

        end
       
        
        function func = checkVariables(func)

            % store function type data
            isSpaceVarying = num2str(Coefficients1d.hasVariable(func,'x'));
            isTimeVarying  = num2str(Coefficients1d.isTimeVarying(func));
            isNonlinear    = num2str(Coefficients1d.isNonlinear(func));
            isCoupled      = num2str(Coefficients1d.isCoupled(func));
            code = strcat(isCoupled, isNonlinear, isTimeVarying, isSpaceVarying);
            code = str2num(code);

            % update function signature based on variable dependence
            switch code

                % no variables
                case 0000
                    func = @(x,t,u,v) func();
                
                % spatially varying only
                case 0001
                    func = @(x,t,u,v) func(x);

                % time varying only
                case 0010
                    func = @(x,t,u,v) func(t);

                % space and time varying
                case 0011
                    func = @(x,t,u,v) func(x,t);

                % nonlinear only
                case 0100
                    func = @(x,t,u,v) func(u);

                % nonlinear and spatially varying
                case 0101
                    func = @(x,t,u,v) func(x,u);

                % nonlinear and time varying
                case 0110
                    func = @(x,t,u,v) func(t,u);

                % nonlinear, space and time varying
                case 0111   
                    func = @(x,t,u,v) func(x,t,u);

                % coupled only
                case 1000
                    func = @(x,t,u,v) func(v);

                % coupled and spatially varying
                case 1001
                    func = @(x,t,u,v) func(x,v);

                % coupled and time varying
                case 1010
                    func = @(x,t,u,v) func(t,v);

                % coupled, space and time varying
                case 1011
                    func = @(x,t,u,v) func(x,t,v);

                % coupled and nonlinear
                case 1100
                    func = @(x,t,u,v) func(u,v);

                % coupled, nonlinear and spatially varying
                case 1101
                    func = @(x,t,u,v) func(x,u,v);

                % coupled, nonlinear and time varying
                case 1110
                    func = @(x,t,u,v) func(t,u,v);

                % coupled, nonlinear, space and time varying
                case 1111
                    func = @(x,t,u,v) func(x,t,u,v);

            end
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
				f = symfun(f);
			end
			
			if isa(f,'sym') || isa(f,'symfun')
				f = matlabFunction(f);
			end
			
			usedVars = Coefficients.getFunctionVariables(f);
			result = any(strcmp(usedVars,checkVar));
		end

    end
end

