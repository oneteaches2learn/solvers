function tx = computeTransmissibility(dx,perm)
%COMPUTETRANSMISSIBILITY computes the transmissibility at each edge

    n    = length(dx);
    temp = zeros(n+1,1);
    temp(1:n) = perm./dx;
    tx   = zeros(n+1,1);

    for i = 2:n
        tx(i) = harmonicMean(temp(i-1),temp(i));
    end
    tx(1) = 2*temp(1);
    tx(n+1) = 2*tx(n);

end



function mean = harmonicMean(x,y)
%HARMONICMEAN(x,y) computes harmonic mean of x and y
    
    mean = 2/(1/x + 1/y);

end


