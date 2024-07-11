% TRUE SOLUTION
u = @(x,t)(sin(pi*x)+0*t);
% SOLVER PARAMETERS
L = 1;
k = 1; c = 0; theta = 1; f = @(x,t)(pi^2*sin(pi*x)+0*t);
por = 1; init = @(x)(sin(pi*x));
% TEST PARAMETERS
n = 4; xc = {}; tGrid = {}; U = {}; V = {}; Norm = {}; Ratio = {}; Log = {};
% X: TEST
for i = 1:n,
nx = 2^(4+(i-1)); nt = 2^(4+2*(i-1)); dx = L/nx; dt = L/nt;
[U{i},xc{i},tGrid{i}] = pennesCCFD1dtime({L,nx},{k,c,theta,f,por},{0,'D',0,'D'},{1/nt,nt,init});
[X,Y] = meshgrid(xc{i},tGrid{i});
V{i} = u(X,Y);
Norm{i} = sqrt(sum(dx*dt*(U{i}-V{i}).^2,'all')); % L2 error
end;
% X: PRINT RESULTS
wait = waitbar(0,'Printing results...');
for i = 1:n-1
Ratio{i} = Norm{i}/Norm{i+1}; Log{i} = log10(Ratio{i})/log10(2);
waitbar(i/n, wait, 'Please wait...');
end;
close(wait);
Norm, Ratio, Log
% T: TEST
for i = 1:n,
nx = 2^(4+2*(i-1)); nt = 2^(4+(i-1)); dx = L/nx; dt = L/nt;
[U{i},xc{i},tGrid{i}] = pennesCCFD1dtime({L,nx},{k,c,theta,f,por},{0,'D',0,'D'},{dt,nt,init});
[X,Y] = meshgrid(xc{i},tGrid{i});
V{i} = u(X,Y);
Norm{i} = sqrt(sum(dx*dt*(U{i}-V{i}).^2,'all')); % L2 error
end;
% T: PRINT RESULTS
wait = waitbar(0,'Printing results...');
for i = 1:n-1
Ratio{i} = Norm{i}/Norm{i+1}; Log{i} = log10(Ratio{i})/log10(2);
waitbar(i/n, wait, 'Please wait...');
end;
close(wait);
Norm, Ratio, Log

surfplotterCCFD(U{4},L,dt,16);
