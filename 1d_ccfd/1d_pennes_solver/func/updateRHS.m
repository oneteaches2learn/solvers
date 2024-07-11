function Q = updateRHS(tx,Q,boundary,dt,t)
%UPDATERHS updates vector F, which represents the right-hand side, based 
%   on the boundary conditions.

    % unpack boundary variables
    leftBC      = boundary{1};
    leftBCtype  = boundary{2};
    rightBC     = boundary{3};
    rightBCtype = boundary{4};

    % store # of intervals
    n = length(tx)-1;

    % apply x = 0 boundary condition
    if leftBCtype == 'D'
		if nargin(leftBC) == 0
        	Q(1) = Q(1) + dt * tx(1) * leftBC();
		elseif nargin(leftBC) == 1
        	Q(1) = Q(1) + dt * tx(1) * leftBC(t);
		end
    elseif leftBCtype == 'N'
		if nargin(leftBC) == 0
			Q(1) = Q(1) - dt*leftBC();
		elseif nargin(leftBC) == 1
			Q(1) = Q(1) - dt*leftBC(t);
		end
    elseif leftBCtype == 'R'
		if nargin(leftBC{2}) == 0
			alpha = leftBC{1};
			sigma = leftBC{2};
			Q(1) = Q(1) + dt*tx(1)/(1 + tx(1)*(1/alpha()))*(sigma());
		elseif nargin(leftBC{2}) == 1
			alpha = leftBC{1};
			sigma = leftBC{2};
			Q(1) = Q(1) + dt*tx(1)/(1 + tx(1)*(1/alpha()))*(sigma(t));
		end
    end
    
    % apply x = L boundary condition
    if rightBCtype == 'D'
		if nargin(rightBC) == 0
			Q(n) = Q(n) + dt * tx(n+1) * rightBC();
		elseif nargin(rightBC) == 1
			Q(n) = Q(n) + dt * tx(n+1) * rightBC(t);
		end
    elseif rightBCtype == 'N'
		if nargin(rightBC) == 0
			Q(n) = Q(n) - dt*rightBC();
		elseif nargin(rightBC) == 1
			Q(n) = Q(n) - dt*rightBC(t);
		end
    elseif rightBCtype == 'R'
		if nargin(rightBC{2}) == 0
			alpha = rightBC{1};
			sigma = rightBC{2};
			Q(n) = Q(n) + dt*tx(n+1)/(1 + tx(n+1)*(1/alpha()))*(sigma());
		elseif nargin(rightBC{2}) == 1
			alpha = rightBC{1};
			sigma = rightBC{2};
			Q(n) = Q(n) + dt*tx(n+1)/(1 + tx(n+1)*(1/alpha()))*(sigma(t));
		end
    end
end
