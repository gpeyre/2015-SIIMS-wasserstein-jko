%%
% Simulation of crowd motion through JKO flow.
% The function to be minimized is
%       f(p) = <w,p> + epsilon*E(p) + h(p)
% w act as an attraction potential, and epsilon is the diffusion strength.

addpath('toolbox/');
addpath('imgaussian/');
set(0,'DefaultFigureWindowStyle','docked');

if not(exist('name'))
    name = 'bump';
    name = 'holes';
    name = 'lena';
    name = 'tworooms';
    name = 'disksquare';
    name = 'hibiscus';
    name = 'roof1';
    name = 'aniso1';
end


porous_mode = 'constant';
porous_mode = 'varying-e';
porous_mode = 'varying-em';
porous_mode = 'varying-m';

%%
% Select either crowd motion model or non-linear (e.g. porous) model.
% This will condition which function h is used.

% grid size
N = 200;
N = 100;

%%
% helpers

normalize = @(x)x/sum(x(:));
% set figure title
setfigname = @(name)set(gcf, 'Name', name, 'NumberTitle','off');
t = linspace(0,1,N); [Y,X] = meshgrid(t,t);
gaussian = @(m,s)exp( -( (X-m(1)).^2+(Y-m(2)).^2 )/(2*s^2) );

%% 
% Load attracting potential and initial density

[w,p0,M,mask] = load_setup(name, N);

%%
% Parameters, here for N=100

[epsilon,tau,gamma,model] = load_default_parameters(name,N);

rep = ['results/' model '/'];
if not(exist(rep))
    mkdir(rep);
end

%%
% Display levelset of w.

w1 = w.*mask; w1 = w1/max(w1(:));
f = cos(150*w1);
figure(1); setfigname('Potential');
clf; 
imageplot(f);
imwrite(rescale(f), [rep name '-potential.png'],'png');

%%
% Display metric if needed.

if not(isempty(M))
    f = cos(80*w1);
    opt.sub = 8;
    opt.color = 'r';
    clf;
    plot_tensor_field(M, f, opt);
    drawnow;
end

%%
% Crowd motion model parameter.

if not(exist('kappa_mult'))
    kappa_mult = 1;
end
kappa = max(p0(:))*kappa_mult; % box constraint on the density

if strcmp(name(1:4), 'roof')==1
    % varying kappa
    i = str2num(name(5));
    kappa = max(p0(:)) * rescale(Y, 1-i/4, 1);
    kappa = ones(N); kappa(:,1:round(2*N/3)) = .1;
    kappa = Inf; % kappa * max(p0(:));
end

%%
% Porous-media model parameter.

if strcmp(model, 'porous') %  not(exist('m_porous'))
    m_porous = 1.4;    
    switch porous_mode
        case 'constant'
            m_field = m_porous*ones(N);
            e_field = 5*ones(N);
        case 'varying-m'
            m_field = rescale(gaussian([.5 .5],.2),1.01,2);
            e_field = ones(N);
            %
            epsilon = .05;
            tau = 10/2;
            gamma = .3;
        case 'varying-e'
            m_porous = 1.1; 
            m_field = m_porous*ones(N);
            e_field = rescale( gaussian([.5 .5],.2), .1, 2);
            %
            epsilon = .05;
            tau = 10/2;
            gamma = .3;
        case 'varying-em'
            m_field = rescale(gaussian([.5 .5],.2),1.01,2);
            e_field = rescale( gaussian([.5 .5],.2), .1, 2);
            %
            epsilon = .05;
            tau = 10/2;
            gamma = .3;
    end
end

%%
% Load the proximal map.

switch  model 
    case 'crowd'
        % crowd motion clamping
        proxh = @(p,sigma)min(p,kappa);  
        options.algorithm = 'projection';
        options.algorithm = 'dykstra';
        paramstr = sprintf('kappa%d',10*kappa_mult);
    case 'porous'
        if isscalar(m_field)
            % constant medium
            ProxPorous = load_porous_prox();
            proxh = @(p,sigma)ProxPorous(p,m_porous,sigma);
        else
            % spacially varying medium
            opts.niter = 20; % newton steps
            proxh = @(p,sigma)porous_prox_newton(sigma*e_field,m_field,p, opts);
        end
        options.algorithm = 'dykstra';
        % proxh = @(p,sigma)p;
        paramstr = porous_mode; % sprintf('m%d',10*m_porous);
end

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
    [blur, Delta, Grad] = blurAnisotropic(M,opt);
    filtIter = 5;
    K = @(x)blur(x,gamma,filtIter);
end

%%
% Test single step.

options.niter = 50; 
options.tol = 1e-15;
options.nsteps = 1;
options.niter_deconvol = 0;
figure(3); setfigname('Iterations');
clf;
[p,Constr] = perform_jko_stepping(K,w,p0, gamma,tau,epsilon,proxh, options);
% display
figure(2); setfigname('Constraints');
clf;
subplot(2,1,1);
plot(1:length(Constr{1}), log10(Constr{1}));
subplot(2,1,2);
plot(1:length(Constr{2}), log10(Constr{2}));
drawnow;

%%
% Perform a deconvolution of the result.

Kdeconv = @(x)blur(x,gamma);
Kdeconv = K;
options.Kdeconv = Kdeconv;
opt.niter = 30;
[q,err] = deconv_richlucy(Kdeconv,p,opt);

figure(4); 
clf; setfigname('Deconvolution');
imageplot({p0 K(p0) p q}, {'p_0' 'K(p_0)' 'p' 'Deconv(p)'}, 2, 2);

%%
% Perform several iterates.

figure(3); setfigname('Iterations');
options.DispFunc = @(p)imageplot(p);
options.DispFunc = @(p)imageplot( format_image(p,mask) ); % 
options.WriteFunc = []; % @(p,i)imwrite(rescale(p), [rep1 num2str(i) '.png'], 'png');
options.niter = 100; 
options.nsteps = 50;
options.tol = 1e-6;
options.niter_deconvol = 0; % 50; % 200;
[p,Constr,plist] = perform_jko_stepping(K,w,p0, gamma,tau,epsilon, proxh, options);

% convert to displayable space/time volume
A0 = reshape(cell2mat(plist), [N N options.nsteps]);
A0(isnan(A0))=0;
Mask = repmat(mask, [1 1 size(A0,3)]);
A = rescale(A0); 
% normalize each frame to max contrast
A(Mask==0) = 0;
A = A - repmat(min(min(A,[],1),[],2), [N N 1]);
A = A ./ repmat(max(max(A,[],1),[],2), [N N 1]);
A(Mask==0) = 1;

% save to video file
write_video(A, [rep name '-' paramstr], 'gif');
% save to folder
nbr = 20;
isvg = round( linspace(1,length(plist),nbr) );
rep1 = [rep name '-' paramstr '/'];
if not(exist(rep1))
    mkdir(rep1);
end
for i=1:length(isvg)
    imwrite( A(:,:,isvg(i)), [rep1 name '-' paramstr '-' num2str(i) '.png'], 'png' );
end

% a = squeeze(sum(sum(A0,1),2)); plot(a); axis tight;

if 0
t = linspace(0,1,N);
[Y,X] = meshgrid(t,t);
d = rescale( atan((Y-.5)*20), .01,.99 );
d = d/d(end/2,end)*p0(end/2,end);

plot([p0(end/2,:); p(end/2,:)]');

imageplot(p-d);
plot(p(end/2,:)-d(end/2,:));
end