%%
% Test for joint evolution of densities with penalization of the sum.

addpath('toolbox/');
addpath('imgaussian/');
set(0,'DefaultFigureWindowStyle','docked');

if not(exist('name'))
    name = 'congestion';
    name = 'tworectangles';
end

rep = ['results/pairwise-sum/'];
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

donormalize = 1;
switch name
    case 'tworectangles'
        s1 = .15; s2 = 1;
        c = [.5 .5];
        p01 = rectange(c, [s1 s2]);
        p02 = rectange(c, [s2 s1]);
        % target flow
        w1 = zeros(N);
        w2 = zeros(N);
        M = [];
        mask = ones(N);
        % coupling is entropy
        coupling = 'entropy';
        kappa_mult = Inf;
        % parameters N=100
        tau = 4;
        gamma = .5;
        % parameters N=200
        tau = 4*2;
        gamma = .5*2;
        addstr = '';
    case {'congestion' 'congestion-sizes'} 
        if strcmp(name, 'congestion')
            shape = { disk disk };
            shape = { square square };
            s = .25*[1 1];
            c = { [.5 .5-s(1)] [.5 .5+s(2)] };
        elseif strcmp(name, 'congestion-sizes')
            shape = { disk disk };
            donormalize = 0;
            s = .25*[1 .7];
            d = .37;
            c = { [.5 .5-s(1)] [.5 .5+s(2)] };
        end
        p01 = shape{1}(c{1}, s(1));
        p02 = shape{2}(c{2}, s(2));
        % target flow 
        w2 = Y; w1 = 1-Y;      
        % coupling is bound
        coupling = 'congestion';
        kappa_mult = 1;
        % parameters N=100
        tau = 200;
        gamma = .4;
        % parameters N=100
        tau = 500;
        gamma = .8;
        addstr = ''; % ['-' num2str(round(10*kappa_mult))];
    otherwise
        error('Unknown setting');
end
vmin = .01;
p01 = p01+vmin;
p02 = p02+vmin;
if donormalize
    p01 = normalize(p01);
    p02 = normalize(p02);
end

rep1 = [rep name addstr '/'];
cmkdir(rep1);


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
% Couplings

switch coupling 
    case 'entropy'
        proxh = @(p)p.^(1/(1+tau/gamma));
    case 'congestion';
        proxh = @(p)min(p,kappa);
    otherwise
        error('Unknown setting');
end

%%
% Run dykstra's algorithm for single stepping.

q1 = p01; q2 = p02; 
%
m1 = exp(-tau/gamma*w1);
m2 = exp(-tau/gamma*w2);
%
multc = @(a,v){a{1}.*v, a{2}.*v};
proxE = @(a)multc(a, proxh(a{1}+a{2})./(a{1}+a{2})  );
prox1 = @(a)proxE({a{1}.*m1, a{2}.*m2});
prox2 = @(b){q1 q2}; 
%
options.DispFunc = {@(p)imageplot(fuse(p{1},p{2})), @(p)[]};
options.niter = 10;
options.Err = {@(p)1, @(p)celldist(p, {q1 q2}) };
options.tol = 1e-8;
figure(2); setfigname('Iterations');
clf;
[a,b,Constr,p_list] = perform_dikstra_scaling(K, [N N 2], { prox1 prox2 }, options);

% next steps
p1 = a{1}.*K(b{1}); 
p2 = a{2}.*K(b{2}); 

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

options.niter = 300;
nsteps = 50; 
options.DispFunc = {@(p)[], @(p)[]};
A1 = []; A2 = []; q1 = p01; q2 = p02;
figure(4); setfigname('Iterations');
for i=1:nsteps
    A1(:,:,i) = q1;
    A2(:,:,i) = q2;
    % clf; imageplot(fuse(im(q1),im(q2))); drawnow;
    clf; imageplot(q1+q2); drawnow;
    % iterations
    prox2 = @(b){q1 q2}; 
    fprintf([num2str(i) ':']);
    [a,b] = perform_dikstra_scaling(K, [N N 2], { prox1 prox2 }, options);
    q1 = a{1}.*K(b{1}); 
    q2 = a{2}.*K(b{2});
end

%%
% Save to files.

% Mask = repmat(mask, [1 1 size(A0,3)]);

A = cat(4,A1,A2,A1*0);
A = permute(A, [1 2 4 3]);
if not(isinf(kappa))
    A = min(A/kappa,1);
else
    A = A - repmat(min(min(A,[],1),[],2), [N N 1 1]);
    rescaleFrame = @(A)A ./ repmat(max(max(A,[],1),[],2), [N N 1 1]);
    A = rescaleFrame(A); A(isnan(A))=0;
end
write_video(A, [rep name addstr], 'gif');
As = rescaleFrame( squeeze(sum(A,3)) ); As(isnan(As))=0;
write_video(As, [rep name addstr '-sum'], 'gif');
write_video(squeeze(A(:,:,1,:)), [rep name addstr '-1'], 'gif');
write_video(squeeze(A(:,:,2,:)), [rep name addstr '-2'], 'gif');


% save to folder
nbr = 20;
isvg = round( linspace(1,nsteps,nbr) );
for i=1:length(isvg)
    a = A(:,:,:,isvg(i));
    a(isnan(a)) = 0;
    imwrite( rescale(a), [rep1 name '-' num2str(i) '.png'], 'png' );
    imwrite( rescale(squeeze(sum(a,3))), [rep1 name '-sum-' num2str(i) '.png'], 'png' );
    imwrite( rescale(squeeze(a(:,:,1))), [rep1 name '-p1-' num2str(i) '.png'], 'png' );
    imwrite( rescale(squeeze(a(:,:,2))), [rep1 name '-p2-' num2str(i) '.png'], 'png' );
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