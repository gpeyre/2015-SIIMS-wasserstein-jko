%% Set up paths

path(path,'mesh_functions/');
path(path,'meshes/');
addpath('toolbox/');

%%
% Helpers

normalize = @(h)h/sum(h(:));

%%
% Read mesh.

name = 'moomoo';
name = 'elephant';

rep = ['results/meshes/' name '/'];
if not(exist(rep))
    mkdir(rep);
end

[X,T] = readOff([name '.off']);
M = getMeshData(X,T,10); % compute 10 LB eigenfunctions for fun
n = M.numVertices;

clf;
showDescriptor(M,M.laplaceBasis(:,5));

%%
% Display heat kernel

s = 1; % vertex # for the dirac location
if strcmp(name, 'elephant')
    s = 15000; % tip of the trompe
end

x = zeros(n,1); x(s)=1;
v = M.vertices(s,:);

tlist = [.01 .05 .1 .4 1 2];
blurSteps = 5;
for i=1:length(tlist)
    t = tlist(i);
    u = blurOnMesh(x,M,t,blurSteps);
    clf; showDescriptor(M,u);
    drawnow;
    saveas(gcf, [rep name '-blur-' num2str(i) '.png'], 'png');
end

%%
% compute a few Gaussians

% farthest point sampling
p = 10; % number of seeded points
t = .1; bs = 5; % blur to comptue far points
slist = [s];
for i=1:p-1
    % compute blur to seed points
    x = zeros(n,1); x(slist) = 1;
    u = blurOnMesh(x,M,t,bs);
    [~,slist(end+1)] = min(u);
end

%%
% Compute input densities

t = .005; bs = 10;
for i=1:p
    x = zeros(n,1); x(slist(i)) = 1;
    Q{i} = blurOnMesh(x,M,t,bs);
end

Qmax = cell2mat(Q); Qmax = max(Qmax(:));
for i=1:p
    Q{i} = min( Q{i}, Qmax*.7);
end

u = sum( cell2mat(Q), 2 );
clf;  showDescriptor(M, u);
saveas(gcf, [rep name '-input.png'], 'png');


%%
% Set up Gaussian blur function for barycenters

% if this gets too small, distances get noisy
blurTime = .001/10; blurSteps = 6;
blur = @(x)blurOnMesh(x,M,blurTime,blurSteps); % faster than pre-factored?
blurTranspose = @(x)blurOnMesh(x,M,blurTime,blurSteps,1);

%%
% Distance.

[d,w0,w1] = convolutionalDistance(Q{1},Q{2},M.areaWeights,blur,blurTranspose);

clf; showDescriptor(M,log(max(w0,1e-10)));

%%
% Barycenters.

options.tol = 1e-6;
options.niter = 100;
options.verb = 2;
options.disp = [];
options.disp = @(x)[]; % showDescriptor(M, x);
w = normalize(ones(p,1));
clf;
[B,u] = wassersteinBarycenter(Q,blur,w, [], options);

clf;  showDescriptor(M, B);
saveas(gcf, [rep name '-isobarycenter.png'], 'png');


%%
% Displacement interpolation.

tlist = [0 .1 .2 .3 linspace(.4,.6,20) .7 .8 .9 1];
for i=1:length(tlist)
    t = tlist(i);
    w = [1-t,t];
    [B,u] = wassersteinBarycenter({Q{1} Q{6}},blur,w, [], options);
    clf;  showDescriptor(M, B);
    drawnow;
    saveas(gcf, [rep name '-inter-' num2str(i) '.png'], 'png');
end
