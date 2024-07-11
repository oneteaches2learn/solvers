function S = constructStiffnessMatrix(dx,tx,C,Por,dt)
%CONSTRUCTSTIFFNESSMATRIX directs the construction of the stiffness 
%   matrix.

    S1 = buildBaseMatrix(tx);
    S2 = spdiags(C'.*dx,0,length(C),length(C));
    S3 = spdiags(Por'.*dx,0,length(Por),length(Por));

    S = dt*S1 + dt*S2 + S3;

end


function A = buildBaseMatrix(tx)
%BUILDBASEMATRIX builds the tridiagonal matrix using transmissibility. 

    n = length(tx)-1;
    A = sparse(n,n);
    
    for i = 2:n
        gl = i-1;
        gr = i;
    
        A(gl,gl) = A(gl,gl) + tx(i);
        A(gl,gr) = A(gl,gr) - tx(i); 
        A(gr,gl) = A(gr,gl) - tx(i);
        A(gr,gr) = A(gr,gr) + tx(i);
    end
end


