clear all

% domain
a = 1;
b = 2;
M = 100;

% coefficients
syms x;
syms t;
syms u;
k     = 1;
c     = x^2 + 0.5 + t;
uStar = x - 2 + t;
por   = x + 1 + t;

% right hand side
f(x,t,u) = sin(pi*x)*exp(t);

% boundary
leftBC    = 0*t;
leftType  = 'D';
rightBC   = 0*t;
rightType = 'D';

% time
dt  = 1;
nt  = 10;
u_0(x) = 0*x;

% assemble inputs
k   = matlabFunction(k);
c   = matlabFunction(c);
uStar = matlabFunction(uStar);
por = matlabFunction(por);
f   = matlabFunction(f);
u_0 = matlabFunction(u_0);
leftBC  = matlabFunction(leftBC);
rightBC = matlabFunction(rightBC);

domain = {a,b,M};
parameters = {k,c,uStar,f,por};
boundary = {leftBC,leftType,rightBC,rightType};
time = {dt,nt,u_0};

[U,xc,tGrid] = pennesCCFD1dTimeNonlinear(domain,parameters,boundary,time);
surfplotterCCFD(U,{xc,tGrid},dt,16); 
