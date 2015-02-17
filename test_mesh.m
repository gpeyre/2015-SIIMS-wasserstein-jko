%%
% Test for JKO on surfaces.

addpath('toolbox/');
addpath('toolbox_mesh/');

rep = 'results/';
if not(exist(rep))
    mkdir(rep);
end

name = 'moomoo';

%%
% Helpers.

normalize = @(h)h/sum(h(:));
setfigname = @(name)set(gcf, 'Name', name, 'NumberTitle','off');

%%
% Mesh loading.

if not(exist('V')) || not(exist('F'))
    nsub = 2; % #subdivision to increase resolution
    [V0,F0] = read_off([name '.off']);
    options.sub_type = 'loop';
    [V,F] = perform_mesh_subdivision(V0, F0, nsub, options);
    %
    M = getMeshData(V',F',1);
    N = M.numVertices;
end

switch name
    case 'moomoo'
        % UP hands: 874, 349
        % MIDLE: 194
        % DOWN hands: 240, 199
        % ARMS: 438 646
        src = [874, 349, 240, 199, 438 646];
        tgt = [292];
        %
        src = [438, 646, 874, 349];
        tgt = [240, 199];
    otherwise
        error('Unknown mesh');
end

if 0
    %% Vertex selection %%
    clf;
    plot_mesh(V0,F0);
    colormap parula(256);
    %
    mesh.vertices = V0';
    idx=selectPoint(gcf, mesh)
end


%%
% Blurring kernel

blurSteps = 5; 
blur = @(u,mu)blurOnMesh(u,M,mu,blurSteps);

%%
% Compute initial density

s = .001*2; % width of the input
U = [];
for i=1:length(src)
    u = blur(dirac(N,src(i)),s);
    U(:,end+1) = rescale( min(u,max(u)*.7) );
end
p0 = normalize( sum(U,2) );

figure(1); setfigname('Initial density');
clf;
myplot_mesh(V,F,p0);
colormap parula(256);

%%
% Compute potential via distance function

options.nb_iter_max = Inf; 
d = perform_fast_marching_mesh(V, F, tgt, options);
w = d/max(d); 

figure(2); setfigname('Potential');
clf;
myplot_mesh(V,F, cos(120*w));
saveas(gcf, [rep name '-potential.png'], 'png');

%%
% 

switch nsub
    case 1
        epsilon = .01/5; % diffusion strength
        tau = .001*2;
        gamma = .01/600;
    case 2
        epsilon = .01/5; % diffusion strength
        tau = .001*3.5;
        gamma = .01/400;
    otherwise
        error('Parameter not specified.');
end
        
%
if not(exist('kappa_mult'))
    kappa_mult = 2;
end
kappa = max(p0(:))*kappa_mult;

%%
% Kernel

using_factorization = 1;
if not(using_factorization)
    K = @(x)blur(x,gamma);
else
    structure = prefactorMeshBlur(M,gamma,blurSteps);
    K = @(x)prefactoredBlur(x,structure,0);
end

figure(3); setfigname('Iterates');
clf;
options.DispFunc = @(p)myplot_mesh(V,F,p); % 
options.WriteFunc = [];
options.niter = 40; 
options.tol = 1e-5;
options.nsteps = 50;
[p,Constr,plist] = perform_jko_stepping(K,w,p0, gamma,tau,epsilon,kappa, options);

clf;
A0 = [];
for i=1:length(plist)
    myplot_mesh(V,F,plist{i});
    drawnow;
    a = getframe(gcf);
    a = double(a.cdata);
    % set background color to white
    bg = a(1,1,:);
    I = find(a(:,:,1)==bg(1) & a(:,:,2)==bg(2) & a(:,:,3)==bg(3) );
    o = size(a,1)*size(a,2);
    a([I; I+o; I+2*o])=255;    
    %
    A0(:,:,:,i) = a;
end
m = 120; % cropping
A = A0(:,m+1:end-m,:,:)/255;

%%
% Save to image file
nbr = 20;
paramstr = sprintf('kappa%d',10*kappa_mult);
isvg = round( linspace(1,length(plist),nbr) );
rep1 = [rep name '-' paramstr '/'];
if not(exist(rep1))
    mkdir(rep1);
end
for i=1:length(isvg)
    imwrite( A(:,:,:,isvg(i)), [rep1 name '-' paramstr '-' num2str(i) '.png'], 'png' );
end

% write_video(rescale(A), [rep name '-' paramstr], 'avi');
write_video(rescale(A), [rep name '-' paramstr], 'gif');

