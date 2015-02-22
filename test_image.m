%%
% Simulation of crowd motion through JKO flow.
% The function to be minimized is
%       f(p) = <w,p> + h(p)
% w act as an attraction potential.

addpath('toolbox/');
addpath('imgaussian/');
addpath('anisotropic/');
set(0,'DefaultFigureWindowStyle','docked');

if not(exist('name'))
    name = 'bump';
    name = 'tworooms-noncvx';
    name = 'aniso3';
    name = 'tworooms';
    name = 'holes';
    name = 'disk1';
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
N = 99;

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

[tau,gamma,model] = load_default_parameters(name,N);

rep = ['results/' model '/'];
if not(exist(rep))
    mkdir(rep);
end

%%
% Display w and metric.

figure(1); setfigname('Potential');
w1 = w.*mask; w1 = w1/max(w1(:));
if not(isempty(M)) && norm(M(:,:,1,2), 'fro')>1e-5
    f = cos(60*w1);
    opt.sub = 8;
    opt.color = 'r';
    clf;
    plot_tensor_field(M, f, opt);
    drawnow;
    saveas(gcf, [rep name '-potential.eps'], 'epsc');
else
    f = cos(150*w1);
    clf;
    imageplot(f);
    imwrite(rescale(f), [rep name '-potential.png'],'png');
end

%%
% Crowd motion model parameter.

if not(exist('kappa_mult'))
    kappa_mult = 1;
end
kappa = max(p0(:))*kappa_mult; % box constraint on the density

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
            tau = 10/2;
            gamma = .3;
        case 'varying-e'
            m_porous = 1.1; 
            m_field = m_porous*ones(N);
            e_field = rescale( gaussian([.5 .5],.2), .1, 2);
            %
            tau = 10/2;
            gamma = .3;
        case 'varying-em'
            m_field = rescale(gaussian([.5 .5],.2),1.01,2);
            e_field = rescale( gaussian([.5 .5],.2), .1, 2);
            %
            tau = 10/2;
            gamma = .3;
    end
end

%%
% Load the proximal map.

switch  model 
    case 'crowd'
        % crowd motion clamping
        proxf = @(p,sigma)min(p .* exp(-sigma*w),kappa);  
        if strcmp(name, 'tworooms-noncvx')
            % non convex clamping
            a = .0001*max(p0(:));
            b = 1*max(p0(:));
            E = @(x)x.*(log(max(x,1e-10))-1);
            Emean = @(a,b)exp( (E(a)-E(b))/(a-b) );
            pthresh = Emean(a,b);
            proxf1 = @(p,sigma) (p<=pthresh) * a + (p>pthresh) * b  ;
            proxf = @(p,sigma)proxf1(p .* exp(-sigma*w),sigma);
            p0 = proxf1(p0,1);
            options.prox_postprocess = 1;  
        end
        paramstr = sprintf('kappa%d',10*kappa_mult);
    case 'porous'
        if isscalar(m_field)
            % constant medium
            ProxPorous = load_porous_prox();
            proxf = @(p,sigma)ProxPorous(p,m_porous,sigma);
        else
            % spacially varying medium
            opts.niter = 20; % newton steps
            proxf = @(p,sigma)porous_prox_newton(sigma*e_field,m_field,p, opts);
        end
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
    opt.laplacian_type = 'fd'; % finite differences
    if isempty(find(mask==0))
        opt.laplacian_type = 'superbases'; % J.M. Mirebeau's method
    end
    [blur, Delta, Grad] = blurAnisotropic(M,opt);
    filtIter = 5;
    K = @(x)blur(x,gamma,filtIter);
end

%%
% Test single step.

options.niter = 10; 
options.tol = 1e-15;
options.nsteps = 1;
figure(3); setfigname('Iterations');
clf;
[p,Constr] = perform_jko_stepping(K,w,p0, gamma,tau,proxf, options);
% display
figure(2); setfigname('Constraints');
clf;
subplot(2,1,1);
plot(1:length(Constr{1}), log10(Constr{1}));
subplot(2,1,2);
plot(1:length(Constr{2}), log10(Constr{2}));
drawnow;

%%
% Perform several iterates.

figure(3); setfigname('Iterations');
options.DispFunc = @(p)imageplot( format_image(p,mask) ); % 
options.WriteFunc = []; % @(p,i)imwrite(rescale(p), [rep1 num2str(i) '.png'], 'png');
options.niter = 100; 
options.nsteps = 50;
options.tol = 1e-6;
[p,Constr,plist] = perform_jko_stepping(K,w,p0, gamma,tau, proxf, options);

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