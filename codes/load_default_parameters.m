function [tau,gamma,model] = load_default_parameters(name,N)

% load_default_parameters - load some default parameters for an experiments
%
%   [tau,gamma,model] = load_default_parameters(name,N);
%
%   Copyright (c) 2014 Gabriel Peyre

model = 'crowd';
switch name
    case {'tworooms' 'tworooms-noncvx'}   % N=100
        tau = 30*3;
        gamma = .4;
        if strcmp(name, 'tworooms-noncvx')
            tau = 70;
        end
    case 'bump'    
        % for N=200
        tau = 20;
        gamma = .3;
    case 'holes' 
        % N=100
        tau = 40;
        gamma = .15;
    case {'disk' 'disk1'}    
        tau = 35;
        gamma = .2;  
    case {'lena' 'hibiscus' 'rand'}
        model = 'porous';
        % N=200
        % m=1.05 ok
        tau = 10/2;
        gamma = .3; 
    case {'aniso0' 'aniso1' 'aniso2' 'aniso3' 'aniso4'}
        % N=100       
        i = str2num(name(6))+1; 
        tau_list = [5 8 15 25 25];
        gamma_list = [.05 .05 .05 .05 .2];
        tau = tau_list(i);
        gamma = gamma_list(i);        
    otherwise        
        error('Unknown setup');
end

end