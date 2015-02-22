function [Delta,Grad] = compute_anisotropic_laplacian(M, diff_type)

% compute_anisotropic_laplacian - compute anisotropic Laplacian using finite differences
%
%   [Delta,Grad] = compute_anisotropic_laplacian(M);
%
%   Delta=Grad'*Diag(M)*Grad is symmetric *positive*.
%
% Copyright (c) 2015 Gabriel Peyre

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


end