function [x,dx,xc] = constructDomain(domain)
%CONSTRUCTDOMAIN(domain) directs the construction of the domain

    x  = setGridPoints(domain);
    dx = setGridLengths(x);
    xc = setGridMidpoints(x,dx);

end


function x = setGridPoints(domain)
%SETGRIDPOINTS(domain) builds a vector of grid points by parsing input type

    % domain input is a cell array, i.e. {L,n} or {a,b,n}, where:
    %   L = right endpoint of domain, n = number of cells
    if isa(domain, "cell")
        if length(domain) == 2
            L = domain{1};
            n = domain{2};
            x = linspace(0,L,n+1);
        elseif length(domain) == 3
            a = domain{1};
            b = domain{2};
            n = domain{3};
            x = linspace(a,b,n+1);
        else
            error 'Incorrect domain input. Cell array should be format {L, n} or {a, b, n}.'
        end
        
    % domain input is a vector
    elseif isa(domain, 'double')
        if length(domain) == 1
            error 'Incorrect domain input. Vector input should have length at least two.'
        else
            x = domain;
        end
    end
end


function dx = setGridLengths(x)
%SETGRIDLENGTHS builds a vector of lengths of the cells of input x

    dx = diff(x);

end


function xc = setGridMidpoints(x,dx)
%SETGRIDMIDPOINTS(x) builds a vector of cell midpoints
    
    n  = length(x)-1;
    xc = x(1:n) + dx/2;

end


