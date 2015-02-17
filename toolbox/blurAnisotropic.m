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
%   Copyright (c) 2014 Gabriel Peyre

options.null = 0;
diff_type = getoptions(options, 'diff_type', 'fwd');
CholFactor = getoptions(options, 'CholFactor', 0);


N = size(M,1);
id = speye(N);

e = ones(N,1);
switch diff_type
    case 'sym'
        % symmetric
        D = spdiags([e -e]/2, [-1 1], N, N);
        D([1 end],:) = 0; % Neumann BC
    case 'fwd'
        % forward
        D = spdiags([e -e], [0 1], N, N);
        D(end,:) = 0; % Neumann BC
    otherwise
        error('Unknown finite differenciation mode.');
end

Dx = kron(D,id);
Dy = kron(id,D);
Grad = [Dx;Dy];

% metric along the diagonal
diagH = @(h)spdiags(h(:), 0, N^2,N^2);
dH = [diagH(M(:,:,2,2)), diagH(M(:,:,1,2)); ...
      diagH(M(:,:,2,1)), diagH(M(:,:,1,1))];
% Laplacian
Delta = Grad'*dH*Grad;

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

