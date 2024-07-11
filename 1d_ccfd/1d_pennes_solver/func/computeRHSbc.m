function Q_bc = computeRHSbc(tx,boundary,t,u,v)
%UPDATERHS updates vector F, which represents the right-hand side, based 
%   on the boundary conditions.

    % unpack boundary variables
    leftBC      = boundary{1};
    leftBCtype  = boundary{2};
    rightBC     = boundary{3};
    rightBCtype = boundary{4};

    % store # of intervals
    n = length(tx)-1;
	Q_bc = sparse(n,1);

    % apply x = 0 boundary condition
    if leftBCtype == 'D'
		if nargin(leftBC) == 0
        	Q_bc(1) = tx(1) * leftBC();
		elseif nargin(leftBC) == 1
        	Q_bc(1) = tx(1) * leftBC(t);
		elseif nargin(leftBC) == 3
        	Q_bc(1) = tx(1) * leftBC(t,u(1),v);
		end
    elseif leftBCtype == 'N'
		if nargin(leftBC) == 0
			Q_bc(1) = -leftBC();
		elseif nargin(leftBC) == 1
			Q_bc(1) = -leftBC(t);
		elseif nargin(leftBC) == 3
			Q_bc(1) = -leftBC(t,u(1),v);
		end
    elseif leftBCtype == 'R'
		if nargin(leftBC{2}) == 0
			alpha = leftBC{1};
			sigma = leftBC{2};
			Q_bc(1) = tx(1)/(1 + tx(1)*(1/alpha(t)))*(sigma());
		elseif nargin(leftBC{2}) == 1
			alpha = leftBC{1};
			sigma = leftBC{2};
			Q_bc(1) = tx(1)/(1 + tx(1)*(1/alpha(t)))*(sigma(t));
		elseif nargin(leftBC{2}) == 3
			alpha = leftBC{1};
			sigma = leftBC{2};
			Q_bc(1) = tx(1)/(1 + tx(1)*(1/alpha(t,u(1),v)))*(sigma(t,u(1),v));
		end
    end
    
    % apply x = L boundary condition
    if rightBCtype == 'D'
		if nargin(rightBC) == 0
			Q_bc(n) = tx(n+1) * rightBC();
		elseif nargin(rightBC) == 1
			Q_bc(n) = tx(n+1) * rightBC(t);
		elseif nargin(rightBC) == 3
			Q_bc(n) = tx(n+1) * rightBC(t,u(n),v);
		end
    elseif rightBCtype == 'N'
		if nargin(rightBC) == 0
			Q_bc(n) = -rightBC();
		elseif nargin(rightBC) == 1
			Q_bc(n) = -rightBC(t);
		elseif nargin(rightBC) == 3
			Q_bc(n) = -rightBC(t,u(n),v);
		end
    elseif rightBCtype == 'R'
		if nargin(rightBC{2}) == 0
			alpha = rightBC{1};
			sigma = rightBC{2};
			Q_bc(n) = tx(n+1)/(1 + tx(n+1)*(1/alpha(t)))*(sigma());
		elseif nargin(rightBC{2}) == 1
			alpha = rightBC{1};
			sigma = rightBC{2};
			Q_bc(n) = tx(n+1)/(1 + tx(n+1)*(1/alpha(t)))*(sigma(t));
		elseif nargin(rightBC{2}) == 3
			alpha = rightBC{1};
			sigma = rightBC{2};
			Q_bc(n) = tx(n+1)/(1 + tx(n+1)*(1/alpha(t,u(n),v)))*(sigma(t,u(n),v));
		end
    end
end
