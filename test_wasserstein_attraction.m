%%
% Test for Wasserstein attraction.

addpath('toolbox/');
addpath('imgaussian/');
set(0,'DefaultFigureWindowStyle','docked');

if not(exist('name'))
    name = 'bump';
    name = 'aniso3';
    name = 'tworooms';
    name = 'holes';
    name = 'disk';
end

rep = ['results/wasserstein-attraction/'];
rep1 = [rep name '/'];
if not(exist(rep))
    mkdir(rep);
end
if not(exist(rep1))
    mkdir(rep1);
end

% grid size
N = 100;

% helpers
celldist = @(C,D)sum( cell2mat(cellfun(@(x,y)norm(x(:)-y(:)), C, D, 'UniformOutput', false)) );
normalize = @(x)x/sum(x(:));
setfigname = @(name)set(gcf, 'Name', name, 'NumberTitle','off');
cellmin = @(A,B)cellfun(@min, A,B);
cellminS = @(kappa,C){min(C{1},kappa), min(C{2},kappa)};

% shapes
t = linspace(0,1,N); [Y,X] = meshgrid(t,t);
gaussian = @(m,s)exp( -( (X-m(1)).^2+(Y-m(2)).^2 )/(2*s^2) );
rectange = @(c,w)double( max( abs(X-c(1))/w(1), abs(Y-c(2))/w(2) )<1 );
square = @(c,s)rectange(c,[s s]);
diamond = @(c,s)double( abs(X-c(1))/s + abs(Y-c(2))/s <1 );
disk = @(m,s) (X-m(1)).^2+(Y-m(2)).^2 <=s^2;
anulus = @(c,s) disk(c,s) - disk(c,s/2);
shapes = {square diamond anulus disk};

%% 
% Load attracting potential and initial density

[w,p0,M,mask] = load_setup(name, N);

%%
% Target distribution r.

tau = .5;
gamma = 0.15;

tau = .7;
gamma = 0.1;
switch name
    case 'tworooms'
        s = .05;
        c = {[.2 .1] [.8 .1] [.8 .9]};
        a = [1 1 .5];
        r = zeros(N);
        for i=1:length(c)
            r = r + a(i)*gaussian(c{i},s);
        end
    case 'holes'
        r = rectange([.5 .1],[.4 .02]);
    case 'disk'
        s = .05;
        c = {[.5 .1] [.9 .5] [.1 .5]};
        a = [1 .5 .5];
        r = zeros(N);
        for i=1:length(c)
            r = r + a(i)*gaussian(c{i},s);
        end
    otherwise
        error('Uknown setup');
end
r = normalize(r+.03);

%%
% Display.

figure(1); setfigname('Distributions');
clf;
imageplot( {format_image(p0,mask) format_image(r,mask)} )

%%
% Crowd motion model parameter.

if not(exist('kappa_mult'))
    kappa_mult = 1;
end
kappa = max(p0(:))*kappa_mult;

%%
% Gibbs kernel

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

%%
% Run dykstra's algorithm for single stepping.

q = p0;
proxs = {@(p)cellminS(kappa,geometric_mean(p,[1 tau])), @(p){q r}};
% proxs = {@(p)p, @(p)p};
options.DispFunc = {@(p)imageplot(format_image(p{1},mask)), @(p)[]};
options.niter = 10;
options.Err = {@(p)celldist(p, {p{end} p{1:end-1}}), @(p)celldist(p,{q r}) };
options.tol = 1e-8;
figure(2); setfigname('Iterations');
clf;
[a,b,Constr,p_list] = perform_dikstra_scaling(K, [N N 2], proxs, options);

p = a{1}.*K(b{1}); % solution
q1 = a{1}.*K(b{1}); r1 = b{2}.*K(a{2});

% final solution
clf;
imageplot({p0 p}, {'p_0' 'p' 'q_1' 'r_1'});

%%
% Errors display.

mydisp = @(x)log10( max(1e-20,x(2:end)/x(2)) );
figure(3); setfigname('Errors');
clf;
subplot(2,1,1);
plot(mydisp(Constr{1})); axis tight;
subplot(2,1,2);
plot(mydisp(Constr{2})); axis tight;

%%
% Run full JKO flow.

options.niter = 200;
nsteps = 50; 
options.DispFunc = {@(p)[], @(p)[]};
A = []; q = p0; 
figure(4); setfigname('Iterations');
for i=1:nsteps
    A(:,:,i) = q;
    clf; imageplot(format_image(q,mask)); drawnow;
    % iterations
    proxs = {@(p)cellminS(kappa,geometric_mean(p,[1 tau])), @(p){q r}};
    [a,b] = perform_dikstra_scaling(K, [N N 2], proxs, options);
    q = a{1}.*K(b{1});    
end

%%
% Save to files.

% convert to displayable space/time volume
A0 = A;
A0(isnan(A0))=0;
Mask = repmat(mask, [1 1 size(A0,3)]);
A = rescale(A0); 
% normalize each frame to max contrast
A(Mask==0) = 0;
A = A - repmat(min(min(A,[],1),[],2), [N N 1]);
A = A ./ repmat(max(max(A,[],1),[],2), [N N 1]);
A(Mask==0) = 1;

% save target distrib
imwrite( rescale(format_image(r,mask)), [rep1 name '-tgt.png'], 'png' );
% save to video file
write_video(A, [rep name], 'gif');
% save to folder
nbr = 20;
isvg = round( linspace(1,nsteps,nbr) );
rep1 = [rep name '/'];
if not(exist(rep1))
    mkdir(rep1);
end
for i=1:length(isvg)
    imwrite( A(:,:,isvg(i)), [rep1 name '-' num2str(i) '.png'], 'png' );
end