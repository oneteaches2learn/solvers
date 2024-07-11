function result = mms_test(uTrue,domain,parameters,boundary,varargin)
	
	% UNPACK PARAMETERS
	% solver parameters
	if length(domain) == 2
		a = 0;
		b = domain{1};
	elseif length(domain) == 3
		a = domain{1}; 
		b = domain{2};
	end
	k     = parameters{1};
	c     = parameters{2};
	uStar = parameters{3};
	leftBC    = boundary{1};
	leftType  = boundary{2};
	rightBC   = boundary{3};
	rightType = boundary{4};

	% MMS test parameters
	if length(varargin) == 0
		base          = 10; 
		startPower    = 1;
		trials        = 4;
		expectedOrder = 2;
		tolerance     = 0.1;
		displayPlot   = 0;
		displayData   = 0;
	else
		base          = varargin{1}{1};
		startPower    = varargin{1}{2};
		trials        = varargin{1}{3};
		expectedOrder = varargin{1}{4};
		tolerance     = varargin{1}{5};
		displayPlot   = varargin{1}{6};
		displayData   = varargin{1}{7};
	end

	% MANUFACTURE SOLUTION
	% manufacture forcing functions
	du = diff(uTrue);
	q  = -k * du;
	f  = diff(q) + c * (uTrue - uStar);

	% convert symbolic functions to function_handles
	u  = matlabFunction(uTrue);
	du = matlabFunction(du);
	k  = matlabFunction(k);
	q  = matlabFunction(q);
	f  = matlabFunction(f);

	% manufacture boundary conditions
	if leftType == 'D'
		leftBC = u(a);
	elseif leftType == 'N'
		leftBC = k() * du(a);
	elseif leftType == 'R'
		alpha   = 1;
		sigma   = u(a) + 1/alpha * q(a);
		leftBC  = [alpha,sigma];
	end
	if rightType == 'D'
		rightBC = u(b);
	elseif rightType == 'N'
		rightBC = -k * du(b);
	elseif rightType == 'R'
		alpha   = 1;
		sigma   = u(b) - 1/alpha * q(b);
		rightBC = [alpha,sigma];
	end
	boundary = {leftBC,leftType,rightBC,rightType};
	parameters = {k,c,uStar};

	% INITIALIZE STORAGE
	xc = {}; U = {}; V = {}; Norm = {}; Ratio = {}; Log = {};

	% TEST
	for i = 1:trials
		j = i - 1;
		x = linspace(a,b,base^(j+startPower) + 1);
		[U{i},xc{i}] = pennesCCFD1d(Domain(x),parameters,boundary,f);
		V{i} = u(xc{i});
		Norm{i} = sqrt(sum(((b-a)/(base^(j + startPower)+1)).*(U{i}-V{i}).^2));
	end

	% RESULTS
	for i = 1:trials-1
		Ratio{i} = Norm{i}/Norm{i+1}; 
		Log{i} = log(Ratio{i}) / log(base); 
	end

	displayData = 0;
	fprintf('hi')
	%if displayData == 1, Norm, Ratio, Log, end
	if displayPlot == 1, plot(xc{trials},U{trials},xc{trials},V{trials}), end

	if (abs(Log{trials-1} - expectedOrder) < tolerance), result = 1; else, result = 0; end

end
