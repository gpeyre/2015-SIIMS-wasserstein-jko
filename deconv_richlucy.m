function [q,err] = deconv_richlucy(K,p,options)

% deconv_richlucy - Richardson-Lucy KL minimization deconvolution
%
%   [q,err] = deconv_richlucy(K,p,options);
%
% Solve for
%   min_q KL(A(q)|p)   where   A=K o diag(1/K(1))
% using a KL-based gradient descent (aka multiplicative descent, 
% aka Richardson-Lucy algorithm).
%
% Copyright (c) 2015 Gabriel Peyre

options.null = 0;
niter = getoptions(options, 'niter', 100);

u = p*0+1;
v = 1 ./ K(u);
% operator to inverse
A  = @(x)K( v.*x );
AS = @(x)v .* K( x );

mysum = @(x)sum(x(:));
KL = @(x,y)mysum( x.*log(x./y) - x + y );
E =@(q)KL(A(q),p);
E0 = E(p);

q = p;
err = [];
for i=1:niter
    q = q .* AS( p ./ A(q) );
    if nargout>1
        err(i) = E(q)/E0;
    end
end

end
