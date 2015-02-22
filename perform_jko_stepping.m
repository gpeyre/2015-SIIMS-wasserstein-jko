function [p,Constr,p_list] = perform_jko_stepping(K,w,q, gamma,tau, proxf , options)

% perform_jko_stepping - perform JKO descent steps
%
%   [p,Constr] = perform_jko_stepping(K,w,q, gamma,tau, proxf, options)
%
%   If options.nsteps=1, solve for
%       JKO_{tau*f}(q) = argmin_p  W_gamma(q,p) + tau*f(p)
%   Otherwise, initialize p=q, and perform options.nsteps iterations of
%       q <- JKO_{tau*f}(q)
%
%   proxf(p,sigma) must implement the KL prox
%       argmin_{p'} KL(p'|p) + sigma*f(p)
%   Set proxf=[] if no function f is used. 
%
%   Here W_gamma is the entropically smooth wasserstein distance, 
%   for a ground metric c, and 
%   where the Gibbs kernel exp(-c/gamma) is implemented by the callback
%   function K(x) (e.g. a Gaussian blurring)
%
%   Copyright (c) 2015 Gabriel Peyre

% iterates are
%    pi = diag(a)*K*diag(b)
% should converge to the constraint
%    (C1) pi^T*1 = b.*K(a) = q  
%    (C2) pi*1 = a.*K(b) = p

niter = getoptions(options, 'niter', 500);
nsteps = getoptions(options, 'nsteps', 1);
verb = getoptions(options, 'verb', 1);
DispFunc = getoptions(options, 'DispFunc', []);
WriteFunc = getoptions(options, 'WriteFunc', []);
tol = getoptions(options, 'tol', 1e-6);
prox_postprocess = getoptions(options, 'prox_postprocess', 0);

if isempty(proxf)
    % nothing 
    proxf = @(p)p;
end

% helpers
mynorm = @(x)norm(x(:));

if nsteps>1
    options.nsteps = 1;
    options.verb = 0;
    p = q;
    p_list = {p};
    for it=1:nsteps-1
        if verb==1
            progressbar(it,nsteps-1); 
        end
        [p,Constr] = perform_jko_stepping(K,w,p, gamma,tau, proxf, options);
        if prox_postprocess==1
            p = proxf(p,1);
        end
        p_list{it+1} = p;
        if not(isempty(DispFunc))
            DispFunc(p);
            drawnow;
        end
        if not(isempty(WriteFunc))
            WriteFunc(p,it);
        end
    end
    return;
end

uu = w*0+1; % vector of ones
%%% projection init %%%
a = uu; b = uu;
%%% Dykstra init %%%
u = uu; u1 = uu; 
v = uu; v1 = uu; 

Constr  = {[] []};
p_list = {};
for i=1:niter
    if verb==1
        progressbar(i,niter);
    end
    for it=1:2
        % save previous iterates
        a1 = a;  b1 = b; % (ell-1) iterates
        u2 = u1; u1 = u; % (ell-2) and  (ell-1) iterates
        v2 = v1; v1 = v; % (ell-2) and  (ell-1) iterates
        if it==2
            a = a1 .* u2;
            b = q ./ K(a);
            Constr{1}(end+1) = mynorm( a.*K(b)-p )/mynorm(q);
        else
            ta1 = a1 .* u2;  
            b = b1 .* v2;
            p = proxf(ta1.*K(b), tau/gamma);
            a = p ./ K(b);            
            Constr{2}(end+1) = mynorm( b.*K(a)-q )/mynorm(q);
        end
        u = u2 .* a1 ./ a;
        v = v2 .* b1 ./ b;
    end
    if nargout>2
        p_list{end+1} = p;
    end
	% tolerance check
    if Constr{1}(end)<tol && Constr{2}(end)<tol
        if verb==1
            progressbar(niter,niter);
        end
        return;
    end
end

end