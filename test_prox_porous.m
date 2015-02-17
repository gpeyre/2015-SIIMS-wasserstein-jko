%%
% Test for the prox of a power function.

rep = 'results/prox/';
if not(exist(rep))
    mkdir(rep);
end
setfigname = @(name)set(gcf, 'Name', name, 'NumberTitle','off');



%%
% Define the prox operator of the associated entropy

[Prox,E] = load_porous_prox();

lw = 2; % linewidth
fs = 25; % font size
qmax = 2;
qlist = linspace(1e-5, qmax,1024)';
mlist = [.5 1 2 5 10];

%%
% Display the gengeralized energy Em

v = []; lgd = {};
for i=1:length(mlist)
    m = mlist(i);
	v(:,end+1) = E(qlist,m);
    lgd{end+1} = ['m=' num2str(mlist(i))];
end

figure(1); setfigname('Entropies');
clf;
plot(qlist, v, 'LineWidth', lw);
axis([0 max(qlist) min(v(:)) .5]);
legend(lgd);
set(gca, 'FontSize', fs);
saveas(gcf, [rep 'entropies.eps'], 'epsc');


%%
% Prox of porous media

figure(2); clf; setfigname('Proxes');

sigma_list = [1 2 3];
ymax = 1.5;
for k=1:length(sigma_list)
    sigma = ones(length(mlist),1)*sigma_list(k);
    v = []; lgd = {};
    for i=1:length(mlist)
        v(:,end+1) = Prox(qlist,mlist(i),sigma(i));
        lgd{end+1} = ['m=' num2str(mlist(i))];
    end
    clf;
    h = plot(qlist, v, 'LineWidth', lw);
    axis([0 max(qlist) 0 ymax]);
    set(gca, 'FontSize', fs);
    legend(lgd, 'Location', 'southeast');
    saveas(gcf, [rep 'prox_sigma' num2str(sigma(i)) '.eps'], 'epsc');
end

return;

%%

myProx = @(s,m,sigma)exp(-sigma*m/(m-1)) * s .* exp( -sigma*m*s.^(m-1));
m = 3;
sigma = 1;
clf; hold on;
plot(qlist,Prox(qlist,m,sigma));
plot(qlist,myProx(qlist,m,sigma), 'r');

