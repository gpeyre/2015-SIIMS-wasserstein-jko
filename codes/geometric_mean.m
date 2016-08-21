function q = geometric_mean(p,lambda)

% prox_barycenters - geometric mean of a set of inputs in a cell array
%
%   q = geometric_mean(p,lambda)
%
%   log(q{m}) = sum_i lambda(i) * log(p{i});
%
%   this is the prox of the barycenter functional
%
%   Copyright (c) 2015 Gabriel Peyre

M = length(p);
lambda = lambda/sum(lambda);

h = ones(size(p{1}));
for m=1:M
    h = h + lambda(m) * log(p{m});
end
h = exp(h);
for m=1:M
    q{m} = h;
end
