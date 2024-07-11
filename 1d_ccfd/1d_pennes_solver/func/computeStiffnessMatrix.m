function M_k = computeStiffnessMatrix(tx)
%CONSTRUCTSTIFFNESSMATRIX(TX) builds the tridiagonal matrix using transmissibility. 

	% instantiate storage
    n = length(tx)-1;
    M_k = sparse(n,n);
    
	% assemble global matrix from local matrices
    for i = 2:n
        gl = i-1;
        gr = i;
    
        M_k(gl,gl) = M_k(gl,gl) + tx(i);
        M_k(gl,gr) = M_k(gl,gr) - tx(i); 
        M_k(gr,gl) = M_k(gr,gl) - tx(i);
        M_k(gr,gr) = M_k(gr,gr) + tx(i);
    end
end

