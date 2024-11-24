clear; close all;
%%
%test: AsgnProb
distr_para.dim  = 1;
distr_para.type = {'uniform'};
distr_para.data = {[0,1]};
rps = (0.005:0.01:0.995)';
n_mcs = 1e7;
asgn_prob = AsgnProb(rps, distr_para, n_mcs);
histogram(asgn_prob,10);
%%
%test: AsgnProb
distr_para.dim  = 2;
distr_para.type = {'uniform','uniform'};
distr_para.data = {[0,1],[0,1]};
x1v = 0.025:0.05:0.975;
x2v = 0.025:0.05:0.975;
[x1m,x2m]=meshgrid(x1v,x2v);
rps = [x1m(:),x2m(:)];
n_mcs = 1e7;
asgn_prob = AsgnProb(rps, distr_para, n_mcs);
histogram(asgn_prob,20);
%%
% test: point selection
distr_para.dim  = 2;
distr_para.type = {'normal','normal'};
distr_para.data = {[0,1],[0,1]};
n_rps           = 200;
[rps,asgn_prob] = GFDiscr_RPS(distr_para,n_rps);

%%
% test: point selection
distr_para.dim  = 2;
distr_para.type = {'normal','normal'};
distr_para.data = {[10,2],[10,2]};
n_rps           = 400;
[rps,asgn_prob] = GFDiscr_RPS(distr_para,n_rps);
%%
% test: point selection
distr_para.dim  = 2;
distr_para.type = {'normal','uniform'};
distr_para.data = {[10,2],[0.03,0.05]};
n_rps           = 200;
[rps,asgn_prob] = GFDiscr_RPS(distr_para,n_rps);
%%
% test: point selection
distr_para.dim  = 2;
distr_para.type = {'uniform','lognormal'};
distr_para.data = {[0.03 0.05],[23.8743,0.3087]};
n_rps           = 1000;
[rps,asgn_prob] = GFDiscr_RPS(distr_para,n_rps);
histogram(asgn_prob,20);
plot(rps(:,1),rps(:,2),'o','LineWidth',1,'MarkerSize',5);
xlabel('\itX_{\rm1}'); xlabel('\itX_{\rm1}');
set(gca, 'FontName','Arial','FontSize',12);

%% REFERENCES
%----------------------------------------------------------------------
%Chen JB, Yang JY, Li J. (2016) A GF-discrepancy for point selection in stochastic seismic response analysis of structures with uncertain parameters. Structural Safety 59: 20-31. https://doi.org/10.1016/j.strusafe.2015.11.001
%Chen JB, Yang JS, Jensen HA. (2020) Structural optimization considering dynamic reliability constraints via probability density evolution method and change of probability measure. Struct Multidisc Optim 62: 2499–2516. https://doi.org/10.1007/s00158-020-02621-4
%Yang JS, Jensen HA, Chen JB. (2022) Structural optimization under dynamic reliability constraints utilizing probability density evolution method and metamodels in augmented input space. Struct Multidisc Optim 65: 107. https://doi.org/10.1007/s00158-022-03188-y
%Yang JS, Chen JB, Jensen HA. (2022) Structural design optimization under dynamic reliability constraints based on the probability density evolution method and highly-efficient sensitivity analysis. Probab Eng Mech 68: 103205. https://doi.org/10.1016/j.probengmech.2022.103205
