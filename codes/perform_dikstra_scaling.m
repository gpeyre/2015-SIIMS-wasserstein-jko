function [a,b,Constr,p_list] = perform_dikstra_scaling(K, dims, proxs, options)

% perform_dikstra_scaling - perform Dykstra's algorithm with scaling functiionals
%
%   [a,b, Constr, p_list] = perform_dikstra_scaling(K, proxs, options)
%
%   Solve
%       min_pi KL(pi,K) + f_1(pi*1) + f_2(pi' * 1)
%   where pi = {pi_m}_{m=1}^M is a collection of couplings
%   so that the computed iterates has the form
%       pi_m = diag(a{m}) * K * diag(b{m})
%
%   proxs{i} for i=1,2 implements prox_{f_i}^KL
%
%   s = size(a{m})=size(b{m}) is the size of the problems.
%   M is the number of couplings
%
%   Copyright (c) 2015 Gabriel Peyre

mynorm = @(x)norm(x(:));
celldist = @(C,D)sum( cell2mat(cellfun(@(x,y)norm(x(:)-y(:)), C, D, 'UniformOutput', false)) );

options.null = 0;
niter = getoptions(options, 'niter', 500);
verb = getoptions(options, 'verb', 1);
tol = getoptions(options, 'tol', 1e-6);
Err = getoptions(options, 'Err', {@(p)Inf, @(p)Inf});
DispFunc = getoptions(options, 'DispFunc', {@(p)0, @(p)0});

%%
% Sizes.

s = dims(1:end-1);
M = dims(end);

%%
% initialization.

U = ones(s); % vector of ones
for m=1:M
    a{m} = U;
end
b = a;
u = a; u1 = a; 
v = a; v1 = a; 

Constr  = {[] []};
p_list = {{} {}};
for i=1:niter
    if verb==1
        progressbar(i,niter);
    end
    for it=1:2 % Left/right marginals
        % save previous iterates
        a1 = a;  b1 = b; % (ell-1) iterates
        u2 = u1; u1 = u; % (ell-2) and  (ell-1) iterates
        v2 = v1; v1 = v; % (ell-2) and  (ell-1) iterates
        % compute tilde a and tilde b
        for m=1:M
            ta1{m} = a1{m} .* u2{m};
            tb1{m} = b1{m} .* v2{m};         
            if it==1
                p1{m} = ta1{m} .* K(tb1{m});
                if nargout>2
                    P1{m} = a1{m} .* K(b1{m}); % true marginals
                end
            else
                p1{m} = tb1{m} .* K(ta1{m});  
                if nargout>2
                    P1{m} = b1{m} .* K(a1{m}); % true marginals   
                end
            end
        end
        % record errors
        if nargout>2
            Constr{it}(end+1) = Err{it}(P1);
        end
        % apply prox
        p = proxs{it}(p1);
        % update of a,b,u,v
        for m=1:M
            if it==1
                a{m} = ta1{m} .* p{m} ./ p1{m};
                b{m} = tb1{m};
            else
                b{m} = tb1{m} .* p{m} ./ p1{m};
                a{m} = ta1{m};                
            end
            %
            u{m} = u2{m} .* a1{m} ./ a{m};
            v{m} = v2{m} .* b1{m} ./ b{m};
        end
        % display
        DispFunc{it}(p); drawnow;
        % save
        p_list{it}{end+1} = p1;
    end
	% tolerance check
    if nargout>2 && Constr{1}(end)<tol && Constr{2}(end)<tol
        if verb==1
            progressbar(niter,niter);
        end
        return;
    end
end

end