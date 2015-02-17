function  p =format_image(p,mask)

% format_image - rescale image for display
%
%   p =format_image(p,mask);
%
%   Copyright (c) 2015 Gabriel Peyre

p(isnan(p)) = 0;
p = rescale(p);
p(mask==0) = 1;
