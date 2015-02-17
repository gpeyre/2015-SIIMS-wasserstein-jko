function [Prox,E] = load_porous_prox(logmin,logmax,k)

% load_porous_prox - KL prox of a power function
%
%   [Prox,E] = load_porous_prox(logmin=-10,logmax=10,k=2048)
%
%   Prox(q,m,sigma) computes prox_{sigma*E(.,m)}^KL(q) 
%   where E(.,m) is the generalized entropy
%       E(s,m) = s * (log(s)-1)   if m==1
%              = s * (s^(m-1)-m)/(m-1)  otherwise
% where
%   prox_{h}^KL(q) = argmin_p KL(p|q) + h(p)
%
%   Copyright (c) 2015 Gabriel Peyre

% log(p) + sigma*m*p^(m-1) = log(q)
% log(p^{1/sigma}) + m*p^(alpa-1) = log(q^{1/sigma})

% r = p^{1/sigma}
% log(r) + m*p^(sigma*(m-1)) = log(q^{1/sigma})

% prox =  = [log + sigma*m*(.)^{m-1}]^{-1} o log
% prox =  = exp o [Id + sigma*m*exp^{(m-1) .}]^{-1} o log

if nargin<1
    logmin = -30;
end
if nargin<2
    logmax = 30;
end
if nargin<3
    k = 1024*12;
end

% prox for m=~1
u = linspace(logmin,logmax,k)';
invMap = @(v,m,sigma)interp1_safe(u + sigma*m* (exp((m-1)*u)-1)/(m-1), u, v);
Prox_m = @(q,m,sigma)exp( invMap(log(q),m,sigma) );
% prox for m==1
Prox_1 = @(q,m,sigma)q.^(1/(1+sigma));

Prox = @(q,m,sigma)ProxSelector(Prox_1,Prox_m, q,m,sigma);
E = @(s,m)entropy(s,m);

end

%%%%
function v = ProxSelector(Prox_1,Prox_m, q,m,sigma)

if m==1
    v = Prox_1(q,m,sigma);
else
    v = Prox_m(q,m,sigma);    
end

end

%%%%
function y = entropy(s,m)

if m==1
    y = s.*(log(max(s,1e-10))-1);
else
	y = s.*(s.^(m-1)-m)/(m-1);
end

end

%%%%
function y1 = interp1_safe(x, y, x1)

if min(x1(:))<x(1) 
    if length(find(isinf(x1(:))))>0
        warning('Inf values');
    end
    warning('Out of bound on porous_prox');
    x1(x1(:)<x(1)) = x(1);    
end
if max(x1(:))>x(end) 
    if length(find(isinf(x1(:))))>0
        warning('Inf values');
    end
    warning('Out of bound on porous_prox');
    x1(x1(:)>x(end)) = x(end);
end

y1 = interp1(x, y, x1);

end