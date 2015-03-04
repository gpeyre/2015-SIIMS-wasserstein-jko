%%
% Test for Wasserstein barycenter, using the generic Dykstra's code.


addpath('toolbox/');
addpath('imgaussian/');
set(0,'DefaultFigureWindowStyle','docked');

% grid size
N = 100;

% helpers
celldist = @(C,D)sum( cell2mat(cellfun(@(x,y)norm(x(:)-y(:)), C, D, 'UniformOutput', false)) );
cellnorm = @(C)sum( cell2mat(cellfun(@(x)norm(x(:)), C, 'UniformOutput', false)) );
normalize = @(x)x/sum(x(:));
setfigname = @(name)set(gcf, 'Name', name, 'NumberTitle','off');

% shapes
t = linspace(0,1,N); [Y,X] = meshgrid(t,t);
gaussian = @(m,s)exp( -( (X-m(1)).^2+(Y-m(2)).^2 )/(2*s^2) );
rectange = @(c,w)double( max( abs(X-c(1))/w(1), abs(Y-c(2))/w(2) )<1 );
square = @(c,s)rectange(c,[s s]);
diamond = @(c,s)double( abs(X-c(1))/s + abs(Y-c(2))/s <1 );
disk = @(m,s) (X-m(1)).^2+(Y-m(2)).^2 <=s^2;
anulus = @(c,s) disk(c,s) - disk(c,s/2);
shapes = {square diamond anulus disk};

% number of inputs 
M = 4;

%%
% Regularization.

gamma = 10; % N=100
gamma = 4; % N=200

%%
% Isotropic metric for the Gibbs kernel.

blur = @(x,mu)imgaussian(x,mu,mu*50);
K = @(x)blur(x,gamma);

%%
% input densities
p0 = {}; s = .2;
d = .15;
vmin = .05;
for m=1:M
    c = d + (1-2*d)*rand(2,1);
    p0{m} = normalize( shapes{mod(m-1,M)+1}(c,s)+vmin );
end

figure(1); setfigname('Inputs');
clf; imageplot(p0);

%%
% Barycentric weights

lambda = normalize(ones(M,1));

%%
% Run dykstra's algorithm

proxs = {@(p)geometric_mean(p,lambda), @(p)p0};
% proxs = {@(p)p, @(p)p};
options.DispFunc = {@(p)imageplot(p{1}), @(p)[]};
options.niter = 100;
options.Err = {@(p)celldist(p, {p{end} p{1:end-1}})/cellnorm(p0), @(p)celldist(p,p0)/cellnorm(p0) };
options.tol = 1e-8;
figure(2); setfigname('Iterations');
clf;
[a,b,Constr,p_list] = perform_dikstra_scaling(K, [N N M], proxs, options);

p = p_list{1}(end); p = p{1};
q = p_list{2}(end); q = q{1};

%%
% Errors

mydisp = @(x)log10( max(1e-20,x(2:end)/x(2)) );
figure(3); setfigname('Errors');
clf;
subplot(2,1,1);
plot(mydisp(Constr{1})); axis tight;
subplot(2,1,2);
plot(mydisp(Constr{2})); axis tight;
