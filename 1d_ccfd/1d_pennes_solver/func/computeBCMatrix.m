function M_bc = computeBCMatrix(tx,boundary,t,u,v)    
%UPDATESTIFFNESSMATRIX updates stiffness matrix S based on the 
%   boundary conditions.
%
%-----------------------------------------------------------------------------%
% Notes
%	(1) Robin boundary conditions, 'R', implement the condition
%
%		dot(q(x),n(x)) = -alpha (u(x) - sigma(x)), x in boundary
%
%	where q is flux, n is a unit outward normal, and sigma can be interpreted
%	as the temperature of some outside fluid. The syntax should be one of:
%
%	 	leftBC = {@()(alpha),@()(sigma)};
%		leftBC = {@(t)(alpha(t)),@(t)(sigma(t)};
%		or you can mix time varying / constant.
%	
%-----------------------------------------------------------------------------%

	if nargin < 3, t = 0; end

    % unpack boundary variables
    leftBC      = boundary{1};
    leftBCtype  = boundary{2};
    rightBC     = boundary{3};
    rightBCtype = boundary{4};

    % initialize storage
    n = length(tx)-1;
	M_bc = sparse(n,n);

    % apply x = 0 boundary condition
    if leftBCtype == 'D'
        M_bc(1,1) = tx(1);
    elseif leftBCtype == 'N'
        ...
    elseif leftBCtype == 'R'
        M_bc(1,1) = tx(1)/(1 + tx(1)*(1/leftBC{1}(t,u(1),v)));
    end
    
    % apply x = L boundary condition
    if rightBCtype == 'D'
        M_bc(n,n) = tx(n+1);
    elseif rightBCtype == 'N'
        ...
    elseif rightBCtype == 'R'
        M_bc(n,n) = tx(n+1)/(1 + tx(n+1)*(1/rightBC{1}(t,u(n),v)));
    end
end
