function [w,p0,M,mask] = load_setup(name, N)

% load_setup - load potential w and initial density p0
%
%    [w,p0,M,mask] = load_setup(name, N);
%
%   Copyright (c) 2015 Gabriel Peyre

t = linspace(0,1,N); [Y,X] = meshgrid(t,t);
gaussian = @(m,s)exp( -( (X-m(1)).^2+(Y-m(2)).^2 )/(2*s^2) );
rectange = @(c,w)double( max( abs(X-c(1))/w(1), abs(Y-c(2))/w(2) )<1 );
disk = @(m,s) (X-m(1)).^2+(Y-m(2)).^2 <=s^2;

normalize = @(x)x/sum(x(:));
doclamp = .7;
switch name
    case {'aniso0' 'aniso1' 'aniso2' 'aniso3' 'aniso4'}
        w = .01 + Y;
        mask = ones(N,N);
        p0 = zeros(N,N) +  + .001;
        s = .03;
        pt = @(r,theta)[.5 .5] + .5*r*[cos(theta*pi+pi/2) sin(theta*pi+pi/2)];
        for d = [.1 .3 .6]
            p0 = p0 + gaussian(pt(.8,d) , s);
            p0 = p0 + gaussian(pt(.8,-d), s);
        end
        % anisotropic metric
        i = str2num(name(6)); 
        eta_list = [1 .5 .2 .1 .03];
        eta = eta_list(i+1);
        c = [.5 .5];
        X = X - c(1); Y = Y - c(2);
        d = max(1e-10, sqrt(X.^2+Y.^2));
        e1 = cat(3, X./d, Y./d);
        e2 = cat(3, -Y./d, X./d);
        M = perform_tensor_recomp(e1,e2,ones(N)*eta,ones(N));
    case {'lena' 'hibiscus' 'rand'}
        if strcmp(name, 'rand')
            rand('seed', 123456);
            p0 = rand(N) + .05;
        else
            p0 = load_image(name, N); p0 = mean(p0,3);
            p0 = rescale(crop(p0,N))+.1;
        end
        w = zeros(N,N);
        mask = ones(N,N);
        M = [];
        doclamp = 0;
    case {'tworooms' 'tworooms-noncvx'}
        doclamp = .8;
        % density
        c = [.3 .8];
        s = .15;
        if strcmp(name, 'tworooms-noncvx')
            s = .12;
        end
        p0 = gaussian(c, s);  
        %p0 = disk(c,.09) + .0001;
        % mask
        s = [.45 .18];
        mask = rectange([.5 .75], s) | ... 
            rectange([.5 .25], s) | ... 
            rectange([.5 .5], [.05 .3]);
        mask = double(mask);    
        % mask
        r = .12;
        mask = abs(Y-.5)>r/2;
        r = .08;
        mask( abs(X-.7)<r/2 ) = 1;
        %
        target = round(N*[.1 .03]);
        
    case {'disk' 'disk1'}
        % density
        p0 = gaussian([.5 .9], .14);        
        % mask
        r = .23;
        mask = double( (X-.5).^2+(Y-.45).^2>=r^2 );
        % target point
        target = round(N*[.5 .03]);
        if strcmp(name, 'disk1')
            w = .5*Y;
        end        
    case 'holes'   
        doclamp = .85;     
        s = .15;
        c = [.3 .8];
        % density
        p0 = gaussian(c, s);
        % mask
        r = .1;
        mask = abs(Y-.5)>r/2;
        h = linspace(0,1,7); h = h(2:end-1);
        r = .05;
        for i=1:length(h)
            mask( abs(X-h(i))<r/2 ) = 1;
        end
        target = round(N*[.9 .03]);
    case {'roof0' 'roof1' 'roof2' 'roof3' 'roof4'}
        p0 = gaussian([.5 .85], .04);  
        p0 = disk([.5 .85], .14)+.01;  
        mask = ones(N);
        w = .01 + Y;
        M = [];
        doclamp = 0;
    case 'bump'
        % potential
        t = linspace(-2,2,N);
        [Y,X] = meshgrid(t,t);
        c = [2 0];
        w = (X-c(1)).^2 + (Y-c(2)).^2 + 5*exp(-5*(X.^2+Y.^2)/2 );
        % density.
        t = linspace(0,1,N)';
        q = 5; % number of band
        r = .7; % witdh of the bands
        a = double( 2*abs(q*mod(t,1/q)-.5)<r );
        p0 = a*a' + .01;
        % mask
        mask = ones(N,N);
        % metric
        M = [];
    otherwise
        error('Unknown setting');
end

f = @(u)normalize( u ); % useful to shift because of masked values
if doclamp>0
    f = @(u)normalize( min(u,max(u(:))*doclamp) );
end
p0 = f(p0.*mask+1e-10);

if not(exist('M'))
    % compute a geodesic metric and geodesic potential
    vmin = 0; 
    M = zeros(N,N,2,2);
    M(:,:,1,1) = mask+vmin; M(:,:,2,2) = mask+vmin;
    % potential: geodesic distance to a target
    if not(exist('target'))
        target = round(N*[.5 .03]);
    end
    if not(exist('w'))
        w = .5*compute_geodesic_heat(target,M);
        w(mask==0) = max(w(:))*5;
    end
end