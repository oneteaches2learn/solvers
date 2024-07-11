function boundary = updateBoundaryCoupledODE(boundary,NameValueArgs)
%UPDATEBOUNDARYCOUPLEDODE(BOUNDARY,NAMEVALUEARGS) updates entries of the
%	BOUNDARY cell array as specified by NAMEVALUEARGS
%
% Author: Tyler Fara					Date: April 4, 2024
%-----------------------------------------------------------------------------%
% Notes
%	(1) After input BOUNDARY, other inputs should be specified as name-value
%	pairs. The syntax is like
%
%		updateBoundaryCoupledODE(boundary,L_type = 'N', L_bc = @()(3.0));
%
%	(2) Updated boundary conditions should be function handles, which can be
%	constant or time varying. The syntax is like
%
%		L_bc = @()(4.0) or L_bc = @(t)(4.0*t)
%
%	(3) To explain the syntax like
%
%		if sum(contains(fieldnames(X),'L_alpha')) ... end
%
%	The variable NameValueArgs is a structure with fields that are named and
%	instantiated when the user inputs something like 'L_alpha=@()(4.0)'. Then
%	NameValueArgs.L_alpha is created, and the function_handle @()(4.0) is
%	stored in that field. Therefore, we only want to attempt to update L_bc if
%	the NameValueArgs.L_alpha field exists, i.e. only if a new L_alpha has been
%	passed to the function. So, for convenience, I first rename NameValueArgs
%	as X. The command fieldnames(X) gives the names, as an array of strings, of
%	the fields of X. The command contains(fieldnames(X),'L_alpha') checks if
%	each string in this array and returns a logical 1 or 0 if that string
%	matches 'L_alpha'. This returns an array of 0's, and possibly a 1 if there
%	is a match. So sum(...) sums over those logical values and returns 1 is
%	'L_alpha' is one of the fields, and 0 if not. 
%
%	(4) It would also be a good idea to check that the user has not attempted
%	to supply L_alpha or L_sigma unless L_type = 'R'. But I'm not going to do
%	that right now, so I'm adding this note for later. 
%-----------------------------------------------------------------------------%

    arguments
		boundary
        NameValueArgs.L_bc
        NameValueArgs.L_type
		NameValueArgs.R_bc
		NameValueArgs.R_type
		NameValueArgs.L_alpha
		NameValueArgs.L_sigma
		NameValueArgs.R_alpha
		NameValueArgs.R_sigma
    end

	% for conveniene, rename NameValueArgs
	X = NameValueArgs;

	% unpack input
	L_bc   = boundary{1};
	L_type = boundary{2};
	R_bc   = boundary{3};
	R_type = boundary{4};
	if L_type == 'R', L_alpha = L_bc{1}; L_sigma = L_bc{2}; end
	if R_type == 'R', R_alpha = R_bc{1}; R_sigma = R_bc{2}; end

	% check inputs
	checkInputs(X,L_type,R_type);

	% update boundary conditions
	if sum(contains(fieldnames(X),'L_alpha')), L_alpha = X.L_alpha; end
	if sum(contains(fieldnames(X),'R_alpha')), R_alpha = X.R_alpha; end
	if sum(contains(fieldnames(X),'L_sigma')), L_sigma = X.L_sigma; end
	if sum(contains(fieldnames(X),'R_sigma')), R_sigma = X.R_sigma; end
	if sum(contains(fieldnames(X),'L_type')), L_type = X.L_type; end
	if sum(contains(fieldnames(X),'R_type')), R_type = X.R_type; end
	if sum(contains(fieldnames(X),'L_bc')), L_bc = X.L_bc; end
	if sum(contains(fieldnames(X),'R_bc')), R_bc = X.R_bc; end

	% repackage updated boundary conditions
	if L_type == 'R', L_bc = {L_alpha,L_sigma}; end
	if R_type == 'R', R_bc = {R_alpha,R_sigma}; end
	boundary = {L_bc,L_type,R_bc,R_type};
end


function checkInputs(X,L_type,R_type)
%CHECKINPUTS(X,L_TYPE,R_TYPE) checks if the user is attempting to change the
%	boundary condition type without also supplying a new boundary condition.
%
%-----------------------------------------------------------------------------%
% notes
%	(1) The logic below is complicated. To explain, line by line:
%
%		1. check if L_type was passed as input
%		2.   if so, check that new L_type is actually different from old L_type
%		3.     if so, check if new L_type is D or N
%	 	4.       if so, check if the user has also supplied a new L_bc
%		5.         if not, throw error 1
%       6.     else, check if new L_type is R
%		7.       if so, check if user has supplied *both* L_alpha and L_sigma
%       8.         if not, throw error 2
%-----------------------------------------------------------------------------%
	err1 = 'You cannot change boundary condition type without also providing a new boundary condition.';
	err2 = 'You cannot change to Robin boundary condition without providing both alpha and sigma.';

	if sum(contains(fieldnames(X),'L_type'))
		if ~strcmp(X.L_type,L_type)
			if (strcmp(X.L_type,'D') || strcmp(X.L_type,'N')) && ...
				~sum(contains(fieldnames(X),'L_bc'))
				error(err1)
			end
			if strcmp(X.L_type,'R') && ...
					~((sum(contains(fieldnames(X),'L_alpha')) + sum(contains(fieldnames(X),'L_sigma')))==2)
				error(err2)
			end
		end
	end

	if sum(contains(fieldnames(X),'R_type'))
		if ~strcmp(X.R_type,R_type)
			if (strcmp(X.R_type,'D') || strcmp(X.R_type,'N')) && ...
				~sum(contains(fieldnames(X),'R_bc'))
				error(err1)
			end
			if strcmp(X.R_type,'R') && ...
					~((sum(contains(fieldnames(X),'R_alpha')) + sum(contains(fieldnames(X),'R_sigma')))==2)
				error(err2)
			end
		end
	end
end


