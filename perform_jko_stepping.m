function [p,Constr,p_list] = perform_jko_stepping(K,w,q, gamma,tau,epsilon, proxh , options)

% perform_jko_stepping - perform JKO descent steps
%
%   [p,Constr] = perform_jko_stepping(K,w,q, gamma,tau,epsilon, proxh, options)
%
%   If options.nsteps=1, solve for
%       JKO_{tau*f}(q) = argmin_p  W_gamma(q,p) + tau*f(p)
%           where f(p) = epsilon*E(p) + <p,w> + h(p)
%   Otherwise, initialize p=q, and perform options.nsteps iterations of
%       q <- JKO_{tau*f}(q)
%
%   proxh(p,sigma) must implement the KL prox
%       argmin_{p'} KL(p'|p) + sigma*h(p)
%   Set proxh=[] if no function h is used (pure diffusion). 
%
%   Here W_gamma is the entropically smooth wasserstein distance, 
%   for a ground metric c, and 
%   where the Gibbs kernel exp(-c/gamma) is implemented by the callback
%   function K(x) (e.g. a Gaussian blurring)
%
%   options.algorithm can be set to:
%       - 'projection': uses iterative Bregman projection.
%           Garanteed to work only for f being an affine constraint, 
%           but in practice seems to work also for bound constraint.
%       - 'dykstra': always converges, but more complicated.
%
%   Copyright (c) 2015 Gabriel Peyre

% iterates are
%    pi = diag(a)*K*diag(b)  and  p
% should converge to the constraint
%    (C1) pi^T*1 = b.*K(a) = q  
%    (C1bis) p<=kappa
%    (C2) pi*1 = a.*K(b) = p

niter = getoptions(options, 'niter', 500);
nsteps = getoptions(options, 'nsteps', 1);
verb = getoptions(options, 'verb', 1);
DispFunc = getoptions(options, 'DispFunc', []);
WriteFunc = getoptions(options, 'WriteFunc', []);
tol = getoptions(options, 'tol', 1e-6);
algorithm = getoptions(options, 'algorithm', 'projection');
niter_deconvol = getoptions(options, 'niter_deconvol', 0);

if isempty(proxh)
    % nothing 
    proxh = @(p)p;
end

tau1 = tau*epsilon;

% compute x^a . y^b
GeomMean = @(x,a, y,b) exp( a.*log(x) + b.*log(y)  );
% helpers
mynorm = @(x)norm(x(:));

if nsteps>1
    options.nsteps = 1;
    p = q;
    p_list = {p};
    for it=1:nsteps-1
        [p,Constr] = perform_jko_stepping(K,w,p, gamma,tau,epsilon, proxh, options);
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
a = uu;
b = uu;
p = exp(-w/epsilon);
%%% Dykstra init %%%
u = uu; u1 = uu; 
v = uu; v1 = uu; 
s = uu; s1 = uu;


alpha = tau1 / (tau1+gamma);
Constr{1}  = []; 
Constr{2} = []; 
Constr{3}  = [];
p_list = {};
for i=1:niter
    if verb==1
        progressbar(i,niter);
    end    
    switch algorithm
        case 'projection'
            % projection on C1
            b = q ./ K(a);
            % p = min(p,kappa);
            p = proxh(p,1/epsilon);
            Constr{1}(end+1) = mynorm( a.*K(b)-p )/mynorm(q);
            % projection on C2
            a = GeomMean( p./K(b),alpha,  a,1-alpha );
            p = GeomMean( p,alpha, a.*K(b),1-alpha );
            Constr{2}(end+1) = mynorm( b.*K(a)-q )/mynorm(q);
            % Constr{3}(end+1) = 0; % max( kappa-p(:) );
        case 'dykstra'
            for it=1:2
                % save previous iterates
                a1 = a; p1 = p; b1 = b; % (ell-1) iterates
                s2 = s1; s1 = s; % (ell-2) and  (ell-1) iterates
                u2 = u1; u1 = u; % (ell-2) and  (ell-1) iterates
                v2 = v1; v1 = v; % (ell-2) and  (ell-1) iterates
                if it==1
                    a = a1 .* u2;
                    b = q ./ K(a);
                    p = proxh(p1 .* s2,1/epsilon);
                    Constr{1}(end+1) = mynorm( a.*K(b)-p )/mynorm(q);
                else
                    ta1 = a1 .* u2;
                    tp1 = p1 .* s2;
                    b = b1 .* v2;
                    a = GeomMean( tp1./K(b),alpha, ta1,1-alpha );
                    p = GeomMean( tp1, alpha, ta1 .* K(b), 1-alpha );
                    Constr{2}(end+1) = mynorm( b.*K(a)-q )/mynorm(q);
                end
                u = u2 .* a1 ./ a;
                v = v2 .* b1 ./ b;
                s = s2 .* p1 ./ p;
            end
        otherwise
            error(['Unknown algorithm ' algorithm]);
    end
	% tolerance check
    if Constr{1}(end)<tol && Constr{2}(end)<tol
        progressbar(niter,niter);
        return;
    end
end

if niter_deconvol>0
    opt.niter = niter_deconvol;
    Kdeconv = getoptions(options, 'Kdeconv', K);
    p_old = p;
    [p,err] = deconv_richlucy(Kdeconv,p,opt);
end

end