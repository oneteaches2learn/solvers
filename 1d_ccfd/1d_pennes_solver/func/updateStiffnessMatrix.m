function S = updateStiffnessMatrix(tx,S,boundary,dt)    

%UPDATESTIFFNESSMATRIX updates stiffness matrix S based on the 
%   boundary conditions.

    % unpack boundary variables
    leftBC      = boundary{1};
    leftBCtype  = boundary{2};
    rightBC     = boundary{3};
    rightBCtype = boundary{4};

    % store # of intervals
    n = length(tx)-1;

    % apply x = 0 boundary condition
    if leftBCtype == 'D'
        S(1,1) = S(1,1) + dt * tx(1);
    elseif leftBCtype == 'N'
        ...
    elseif leftBCtype == 'R'
        S(1,1) = S(1,1) + dt*tx(1)/(1 + tx(1)*(1/leftBC{1}()));
    end
    
    % apply x = L boundary condition
    if rightBCtype == 'D'
        S(n,n)   = S(n,n) + dt * tx(n+1);
    elseif rightBCtype == 'N'
        ...
    elseif rightBCtype == 'R'
        S(n,n) = S(n,n) + dt*tx(n+1)/(1 + tx(n+1)*(1/rightBC{1}()));
    end
end

