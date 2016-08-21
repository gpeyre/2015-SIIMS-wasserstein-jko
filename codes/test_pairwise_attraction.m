%%
% Test for joint evolution of densities with pairwise attraction

addpath('toolbox/');
addpath('imgaussian/');
set(0,'DefaultFigureWindowStyle','docked');

if not(exist('name'))
    name = 'twobumps';
end

rep = ['results/pairwise-attraction/'];
cmkdir(rep);

% grid size
N = 100*2;

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

switch name
    case 'twobumps'
        shape = square;
        shape = gaussian;
        s = .06;
        s = .1;
        c = { [.2 .8] [.8 .8]};
        p01 = shape(c{1}, s);
        p02 = shape(c{2}, s);
        % target flow
        w1 = Y;
        w2 = Y;
        M = [];
        mask = ones(N);
    otherwise
        error('Unknown setting');
end
vmin = .01;
p01 = normalize(p01+vmin); 
p02 = normalize(p02+vmin); 

% N=100 param
tau = 100;
gamma = 2;
% attraction
alpha = .0006; % strong
alpha = .0003;
alpha = .0001;
alpha = .0006; % strong

% N=200 param
tau = 300;
gamma = 3;
tau = 200;
gamma = 5;
% attraction
alpha = .0006; % strong
alpha = .0003;
alpha = .0001;

im = @(p)format_image(p,mask);
fuse = @(p1,p2)cat(3,im(p1),im(p2));

%%
% Display.

figure(1); setfigname('Distributions');
clf;
imageplot( fuse(p01,p02) );

%%
% Crowd motion model parameter.

if not(exist('kappa_mult'))
    kappa_mult = 10;
end
kappa = max([p01(:); p02(:)])*kappa_mult;

%%
% Gibbs kernel

[blur,K] = load_kernel(M,mask,gamma);

%%
% Run dykstra's algorithm for single stepping.

proxh = @(p)min(kappa,p);
q1 = p01; q2 = p02; 
geom = @(p,lambda)exp( lambda(1)*log(p{1}) + lambda(2)*log(p{2}) );
%
scale = @(a,p){a.*p{1}, a.*p{2}};
lambda1 = normalize([1 tau*alpha]);
remap1 = @(u,v){u u v};
lambda2 = normalize([tau*alpha 1]);
remap2 = @(u,v){u v v};
m1 = exp(-tau/gamma*w1 / (1+alpha*tau));
m2 = exp(-tau/gamma*w2 / (1+alpha*tau));
%
prox1 = @(a)remap1( proxh( m1 .* geom( {a{1:2}}, lambda1) ), q2 );  
prox2 = @(b)remap2( q1, proxh( m2 .* geom( {b{2:3}},lambda2) ) ); 
%
options.DispFunc = {@(p)imageplot(im(p{1})), @(p)[]};
options.niter = 10;
options.Err = {@(p)celldist(p, {p{2} p{1} q2}), @(p)celldist(p, {q1 p{2} p{3}}) };
options.tol = 1e-8;
figure(2); setfigname('Iterations');
clf;
[a,b,Constr,p_list] = perform_dikstra_scaling(K, [N N 3], { prox1 prox2 }, options);

% next steps
p1 = a{1}.*K(b{1}); 
p2 = b{3}.*K(a{3}); 

% final solution
clf;
imageplot( {fuse(p01,p02), fuse(p1,p2)} );

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

options.niter = 100;
nsteps = 40; 
options.DispFunc = {@(p)[], @(p)[]};
A1 = []; A2 = []; q1 = p01; q2 = p02;
figure(4); setfigname('Iterations');
for i=1:nsteps
    A1(:,:,i) = q1;
    A2(:,:,i) = q2;
    clf; imageplot(fuse(im(q1),im(q2))); drawnow;
    % iterations
    prox1 = @(a)remap1( proxh( m1 .* geom( {a{1:2}}, lambda1) ), q2 );  
    prox2 = @(b)remap2( q1, proxh( m2 .* geom( {b{2:3}},lambda2) ) ); 
    fprintf([num2str(i) ':']);
    [a,b] = perform_dikstra_scaling(K, [N N 3], { prox1 prox2 }, options);
    q1 = a{1}.*K(b{1}); 
    q2 = b{3}.*K(a{3}); 
end

%%
% Save to files.

% Mask = repmat(mask, [1 1 size(A0,3)]);

addstr = num2str(round(alpha*100000));

A = cat(4,A1,A2,A1*0);
A = permute(A, [1 2 4 3]);
A = A - repmat(min(min(A,[],1),[],2), [N N 1 1]);
A = A ./ repmat(max(max(A,[],1),[],2), [N N 1 1]);
A(isnan(A))=0;
write_video(A, [rep name '-' addstr], 'gif');


rep1 = [rep name '-' addstr '/'];
cmkdir(rep1);

% save to folder
nbr = 20;
isvg = round( linspace(1,nsteps,nbr) );
for i=1:length(isvg)
    imwrite( A(:,:,:,isvg(i)), [rep1 name '-' num2str(i) '.png'], 'png' );
end

return;

% convert to displayable space/time volume
A0 = A;
A0(isnan(A0))=0;
A = rescale(A0); 
% normalize each frame to max contrast
A(Mask==0) = 0;
A = A - repmat(min(min(A,[],1),[],2), [N N 1]);
A = A ./ repmat(max(max(A,[],1),[],2), [N N 1]);
A(Mask==0) = 1;


% save to video file
write_video(A, [rep name], 'gif');