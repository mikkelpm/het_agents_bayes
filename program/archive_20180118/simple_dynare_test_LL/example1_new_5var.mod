// Example 1 from Collard's guide to Dynare
var y, c, k, a, b;
varexo e, u;

parameters beta, rho, alpha, delta, theta, psi, tau, std_e, std_u;

alpha = 0.36;
rho   = 0.95;
tau   = 0.025;
beta  = 0.99;
delta = 0.025;
psi   = 0;
theta = 2.95;
std_e = 0.009;
std_u = 0.009;


model;

k = beta*(((exp(b)*c)/(exp(b(+1))*c(+1)))
    *(exp(b(+1))*alpha*y(+1)+(1-delta)*k));
y = exp(a)*(k(-1)^alpha)*(((1-alpha)*y/c/theta)^((1-alpha)/(1+psi)));
k = exp(b)*(y-c)+(1-delta)*k(-1);
a = rho*a(-1)+tau*b(-1) + std_e*e;
b = tau*a(-1)+rho*b(-1) + std_u*u;
end;

varobs y, c;

initval;
y = 1.08068253095672;
c = 0.80359242014163;
//h = 0.29175631001732;
k = 11.08360443260358;
a = 0;
b = 0;
e = 0;
u = 0;
end;

shocks;
var e = 1;
var u = 1;
end;

stoch_simul(irf=0, order=1);
