function myplot_mesh(V,F,C)

% myplot_mesh - plot wrapper
%
%   myplot_mesh(V,F,C);
%
% Copyright (c) 2015 Gabriel Peyre

clf;
plot_mesh(V,F, struct('face_vertex_color', rescale(C)));
colormap parula(256);