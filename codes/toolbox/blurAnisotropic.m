function [Blur, Delta,Grad] = blurAnisotropic(M, options)

% blurAnisotropic - load anisotropic heat kernel
%
%   [Blur, Delta, Grad] = blurAnisotropicblurAnisotropic(M, options);
%
%   M is an (N,N,2,2) representing the 2D Riemannian manifold.
%
%   Blur(u,t,filtIter) iterated filtIter times
%       u = (Id + t*(Delta))\u
%
%   Set options.CholFactor>0 to use a cholesky pre-factorization of the
%   matrices involved. In this case, options.CholFactor is equal to 
%   the time t involved.
%
%   You can set options.laplacian_type to 'fd' (finite differences)
%   or 'superbases'.
%
%   Copyright (c) 2014 Gabriel Peyre

options.null = 0;
diff_type = getoptions(options, 'diff_type', 'fwd');
CholFactor = getoptions(options, 'CholFactor', 0);
laplacian_type = getoptions(options, 'laplacian_type', 'fd');

N = size(M,1);

switch laplacian_type
    case 'fd'
        [Delta,Grad] = compute_anisotropic_laplacian(M, diff_type);
    case 'superbases'
        Grad = [];
        M1 = cat(3, M(:,:,1,1), M(:,:,1,2), M(:,:,2,2));
        Delta = DiffusionSparseMatrix(M1);
    otherwise
        error('Uknown Laplacian type');
end

    
    
% load blurring kernel
if CholFactor==0
    Blur = @(u,t,filtIter)heat_iter(Delta, t, filtIter, u);
else
	t0 = CholFactor;
    % perform cholesky factorization of A
    Id = speye(N*N,N*N);
    A = Id + t0*Delta;
    I = symamd(A);
    R = chol(A(I,I));
    %
    Blur = @(u,t,filtIter)heat_iter_chol(R,I,t, t0, filtIter, u);    
end

end

%%%
function u = heat_iter(Delta, t, filtIter, u)

s = size(u);
u = u(:);
Id = speye(size(u,1));
for i=1:filtIter
    u = (Id + t*Delta)\u;
end
u = reshape(u,s);

end

%%%
function u = heat_iter_chol(R,I,t, t0, filtIter, u)
if t~=t0
    warning('When using CholFactor>0, it must match t.');
end
s = size(u);
u = u(:);
u = u(I); % re-order
for i=1:filtIter
    u = R\(R'\u);
end
u(I) = u; % re-order
u = reshape(u,s);
end

