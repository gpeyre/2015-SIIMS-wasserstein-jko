function [blur,K] = load_kernel(M,mask,gamma)

% load_kernel - load Gibbs convolution-like kernel
%
%   [blur,K] = load_kernel(M,gamma);
%
% Copyright (c) 2015 Gabriel Peyre

if isempty(M)
    % isoptropic metric
    blur = @(x,mu)imgaussian(x,mu,mu*50);
    K = @(x)blur(x,gamma);
else
    % anisotropic metric
    using_factorization = 1;
    if not(using_factorization)
        opt.CholFactor = 0;
    else
        opt.CholFactor = gamma;        
    end
    opt.laplacian_type = 'fd'; % finite differences
    if isempty(find(mask==0))
        opt.laplacian_type = 'superbases'; % J.M. Mirebeau's method
    end
    [blur, Delta, Grad] = blurAnisotropic(M,opt);
    filtIter = 5;
    K = @(x)blur(x,gamma,filtIter);
end

end