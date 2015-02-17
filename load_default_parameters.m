function [epsilon,tau,gamma,model] = load_default_parameters(name,N)

% load_default_parameters - load some default parameters for an experiments
%
%   [epsilon,tau,gamma,model] = load_default_parameters(name,N);
%
%   Copyright (c) 2014 Gabriel Peyre

model = 'crowd';
switch name
    case 'tworooms'     % N=100
        epsilon = .01/6; % diffusion strength
        tau = 60;
        gamma = 1/2;
        % dystra
        epsilon = .01/2;
        tau = 30;
        gamma = 1/2;
    case 'holes'   
        epsilon = .01/4; % diffusion strength
        tau = 35;
        gamma = 1/2;
    case 'bump'    
        % for N=200
        epsilon = .03; 
        tau = 10;
        gamma = .5;
        % N=100
        epsilon = .01/2; % diffusion strength
        tau = 35;
        gamma = .5/2;
    case {'disk' 'disk1'}    
        epsilon = .01/8; % diffusion strength
        tau = 30;
        gamma = 1/2; 
    case 'disksquare' 
        model = 'porous';
        % m=1
        epsilon = .09; 
        tau = 2;
        gamma = .5;  
        % m=1
        epsilon = .2; 
        tau = 1/8;
        gamma = .5;
        % m=4
        epsilon = .09;
        epsilon = .1; 
        tau = 10*20;
        gamma = 1;   
    case {'lena' 'hibiscus' 'rand'}
        model = 'porous';
        % N=200
        % m=1.05 ok
        epsilon = .04; 
        tau = 10/2;
        gamma = .3;  
        % m=1.05 ok
        epsilon = .04; 
        tau = 10*6;
        gamma = .3;  
        % m=1.4 ok
        epsilon = .04; 
        tau = 10*2;
        gamma = .3;  
        % m=4
        epsilon = .01*5; 
        tau = 10*2;
        gamma = .4;
        % m=4
        epsilon = .03; 
        tau = 10*4;
        gamma = .8;
        % m=8
        epsilon = .03; 
        tau = 10*4;
        gamma = .8;
        % varying_m
        epsilon = .05; 
        tau = 10/4;
        gamma = .3;
    case {'roof0' 'roof1' 'roof2' 'roof3' 'roof4'}    
        epsilon = .01; % diffusion strength
        tau = 30*16;
        gamma = .7; 
        %
        epsilon = .01/4; % diffusion strength
        tau = 30*16/2;
        gamma = .7; 
    case {'aniso0' 'aniso1' 'aniso2' 'aniso3'}
        % N=100
        epsilon = .01/6; % diffusion strength
        tau = 15;
        gamma = .3; 
    otherwise        
        error('Unknown setup');
end

end