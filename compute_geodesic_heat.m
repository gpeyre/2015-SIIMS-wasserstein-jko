function D = compute_geodesic_heat(r,M)

% compute_geodesic_heat - compute geodesic distance
%
%   D = compute_geodesic_heat(r,M);
%
%   D is the geodesic distance to r for the metric M using the "Geodesic in Heat" method.
%
%   Copyright (c) 2015 Gabriel Peyre

N = size(M,1);
[blur, Delta, Grad] = blurAnisotropic(M);

% domain
mask = double( sum(sum(M,3),4)>1e-10 );

filtIter = 5;
sigma = 30;

% dirac
u = zeros(N,N); u(r(1),r(2)) = 1;
% nabla(K*u)
v = blur(u,sigma,filtIter);
v = reshape(Grad*v(:), [N*N 2]);
d = sqrt(sum(v.^2,2));
v = v .* repmat(mask(:) ./ max(d,1e-15), [1 2]);
warning off;
D = -Delta \ ( Grad'*v(:) );
warning on;
D = reshape(D, [N N]);
% shift in [0,1]
D = rescale( (D - D(r(1),r(2))) .* mask );