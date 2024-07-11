clear all

v_0 = 1;
mu = 0.5;
dt = 0.01;
N = 100;
c_B = 1;
v_Underbar = NaN;
v_Overbar = 1;
t = [0:dt:dt*100];
U = sin(pi*t).*ones(10,length(t)) + 0.75;
dx = 0.1 * ones(1,10);

v = zeros(1,N+1);
v(1) = v_0;

for i = 2:N+1
	c_B = i*1;
	S = @(x)(integralAverage(x,dx));
	R = @(x,y)(resolvent(x,dt,c_B,v_Underbar,v_Overbar));
	v(i) = bodyTempODE(mu,dt,U(:,i),v(i-1),R,S);
end

plot(t,v);


C_B = c_b(Viter(1));
R(v_in) = @(v_in)(resolvent(v_in,dt,C_B,v_Underbar,v_Overbar);
S(u) = @(u)(integralAverage(u,dx,weight));
v(u,v) = @(u,v)(bodyTempODE(mu,dt,u,v(i-1),R,S));

Viter(2) = computeFunction(v,Uiter(1,:),Viter(1));
