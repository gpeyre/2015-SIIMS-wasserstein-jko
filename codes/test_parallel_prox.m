%%
% test for parallel resolution of the prox

addpath('toolbox/');

N = 100;

%%
% m should not be too large ...
sigma = rand(N,1)*4;
m = 1 + rand(N,1)*4;
s = rand(N,1);

options.niter = 40;
[psi,err,Err] = porous_prox_newton(sigma,m,s, options);

plot(log10(Err'));
I = find(Err(:,end)>1e-3);