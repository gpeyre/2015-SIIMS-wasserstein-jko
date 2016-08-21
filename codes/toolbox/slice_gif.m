function slice_gif(rep, Q)

% slice_gif - slice animated .gif figures in a directory to output set of images
%
%   slice_gif(rep, Q);
%
%   Copyright (c) 2015 Gabriel Peyre

if nargin<1
    rep = 'results/'; 
end
if nargin<2
    Q = 20; % #images
end
    
% addpath('toolbox/');

a = dir([rep '*.gif']);
for i=1:length(a)
    progressbar(i,length(a));
    str = a(i).name;
    str0 = str(1:end-4);
    rep1 = [rep str0 '/'];
    if not(exist(rep1))
        mkdir(rep1);
    end
    [U,c] = imread([rep str]);
    if norm(std(c'))<1e-3
        % BW
        V = U;
    else
        % color
        V = zeros([size(U,1) size(U,2) 3 size(U,4)]);
        for m=1:3
            V(:,:,m,:) = reshape( c(U+1,m), size(U) );
        end
    end
    P = size(V,4);
    isvg = round( linspace(1,P,Q) );
    for k=1:Q
        u = double(V(:,:,:,isvg(k)));
        imwrite( rescale(u), [rep1 str0 '-' num2str(k) '.png'], 'png' );
    end
end