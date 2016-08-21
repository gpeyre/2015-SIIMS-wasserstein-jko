function [psi,err,Err] = porous_prox_newton(sigma,m,s, options)

% porous_prox_newton - compute the proximal operator of the porous medium using Newton iterations
%
%   [psi,err,Err] = porous_prox_newton(sigma,m,s, options)
%
% Allows to solve for the proximal operator in parallel.
% This means that (m,sigma) can be either scalars or vectors of same size as s.  
%
% computes prox_{sigma*E(.,m)}^KL(s) 
%   where E(.,m) is the generalized entropy
%       E(s,m) = s * (log(s)-1)   if m==1
%              = s * (s^(m-1)-m)/(m-1)  otherwise
% where
%   prox_{h}^KL(q) = argmin_p KL(p|q) + h(p)
%
%   Warning: works only for m>1.
%
%   Copyright (c) 2015 Gabriel Peyre

niter = getoptions(options, 'niter', 50);
tol = getoptions(options, 'tol', 1e-10);

%%
% Function to solve for

f  = @(phi,t,sigma,m) phi-t + m .* sigma .* ( exp((m-1).*phi) - 1 ) ./ (m-1);
df = @(phi,t,sigma,m) 1 + m.*sigma.*exp( (m-1).*phi );

% goes to log domain
t = log(s);
phi = t*0;
err = []; Err = [];
for i=1:niter
    F = f(phi,t,sigma,m);
    dF = df(phi,t,sigma,m);
    err(end+1) = mean(abs(F(:)));
    if nargout>2
        Err(:,end+1) = abs(F(:));
    end
    if err(end)<tol
        break;
    end
    phi = phi - F ./ dF;
end
psi = exp(phi);
if err(end)>1e-6
    warning('Maybe Newton did not converge ...');
end
