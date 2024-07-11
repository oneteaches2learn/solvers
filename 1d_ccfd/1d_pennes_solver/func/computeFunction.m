function F = computeFunction(f,varargin)
%COMPUTEFUNCTION evaluates f(x) and stores result as vector F or evalutes
%   f(x,t) and stores result as matrix F.
    if nargin == 2
        x = varargin{1};
        F = evaluateOneInput(f,x);
    elseif nargin == 3
        x = varargin{1};
        t = varargin{2};
        F = evaluateTwoInputs(f,x,t);
	elseif nargin == 4
		x = varargin{1};
		t = varargin{2};
		u = varargin{3};
		F = evaluateThreeInputs(f,x,t,u);
	elseif nargin == 5
		x = varargin{1};
		t = varargin{2};
		u = varargin{3};
		v = varargin{4};
		F = evaluateFourInputs(f,x,t,u,v);
    else
        error 'Too many arguments passed. Expected two, three, four, or five arguments.'

    end
end


function F = evaluateOneInput(f,x)
%EVALUATEONEINPUT evaluates f at x and stores the result in vector F.
    % parse class of f and evaluate f(x) to obtain F
    if isa(f,"double")
        F = f*ones(1,length(x));
    elseif isa(f,"function_handle")
		if nargin(f) == 0
			F = f() * ones(1,length(x));
		elseif nargin(f) == 1
        	F = f(x);
		else
			error('Incorrect number of inputs')
		end
    elseif isa(f,"symfun")
        F = f(x);
        F = double(F);
    else
        error 'Incorrect function type. Expected type is double, symfun, or function_handle.'
    end
end


function F = evaluateTwoInputs(f,x,t)
%EVALUATETWOINPUTS parses function class, evaluates f(x,t), and stores 
%   result in matrix F.

    % parse class of f and evaluate f(x) to obtain F
    if isa(f,"double")
        F = f*ones(1,length(x));
    elseif isa(f,"function_handle")
		if nargin(f) == 0
			F = f()*ones(1,length(x));
		else
        	F = f(x,t);
		end
    elseif isa(f,"symfun")
        F = f(x,t);
        F = double(F);
    else
        error 'Incorrect function type. Expected type is double, symfun, or function_handle.'
    end
end

function F = evaluateThreeInputs(f,x,t,u)
%EVALUATETHREEINPUTS parses function class, evaluates f(x,t,u), and stores 
%   result in matrix F.

    % parse class of f and evaluate f(x) to obtain F
    if isa(f,"double")
        F = f*ones(1,length(x));
    elseif isa(f,"function_handle")
		if nargin(f) == 0
			F = f()*ones(1,length(x));
		else
        	F = f(x,t,u);
		end
    elseif isa(f,"symfun")
        F = f(x,t,u);
        F = double(F);
    else
        error 'Incorrect function type. Expected type is double, symfun, or function_handle.'
    end
end

function F = evaluateFourInputs(f,x,t,u,v)
%EVALUATETHREEINPUTS parses function class, evaluates f(x,t,u), and stores 
%   result in matrix F.

    % parse class of f and evaluate f(x) to obtain F
    if isa(f,"double")
        F = f*ones(1,length(x));
    elseif isa(f,"function_handle")
		if nargin(f) == 0
			F = f()*ones(1,length(x));
		else
        	F = f(x,t,u,v);
			if length(F) == 1
				F = F * ones(size(x));
			end
		end
    elseif isa(f,"symfun")
        F = f(x,t,u,v);
        F = double(F);
    else
        error 'Incorrect function type. Expected type is double, symfun, or function_handle.'
    end
end

