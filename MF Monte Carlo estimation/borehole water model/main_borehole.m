% computes multifidelity sensitivity index and variance estimates for the
% McCormick function and compares to the high-fidelity Monte Carlo estimate

% PAPERS
% E. Qian, B. Peherstorfer, D. O'Malley, V. Vesselinov, and K. Willcox
% Multifidelity estimation of variance and sensitivity indices
% SIAM/ASA Journal on Uncertainty Quantification, 6(2):683-706, 2018.
%
% G. Cataldo, E. Qian, and J. Auclair
% Multifidelity uncertainty quantification and model validation of
% large-scale multidisciplinary systems
% Journal of Astronomical Telescopes, Instruments and Systems, 2022
%
% AUTHORS
% Elizabeth Qian (elizqian@alum.mit.edu) 
% Giuseppe Cataldo
%
% LAST UPDATED
% 18 May 2022

%% SETUP
clear; close all
%addpath('../mfgsa')

% define high-fidelity and low-fidelity models
fcns{1} = @(Z) model1(Z);   % high-fidelity
fcns{2} = @(Z) model2(Z);   % low-fidelity
fcns{3} = @(Z) model3(Z);   % lowest-fidelity

w = [1; 0.2; 0.05];       % assign model weights/costs
vec = [1 1 1];              % says each model is vectorized

d = 8;          % dimension of uncertain input

budget = 4000;   % define computational budget

% if analytical statistics are not available, estimate them. otherwise load
% true statistics from file
estimate = true;   
if estimate
    n_estimate = 100;
    stats = estimate_statistics(fcns,n_estimate);
else
    load('truestats.mat');
end


%% COMPUTE MULTIFIDELITY GLOBAL SENSITIVITY ANALYSIS REPLICATES
n_reps = 100;    % number of replicates 
avg   = zeros(n_reps,2);    vr    = zeros(n_reps,2);
mc_sm = zeros(n_reps,d);    mc_st = zeros(n_reps,d);
mf_sm = zeros(n_reps,d);    mf_st = zeros(n_reps,d);

method = 'Owen';

for n = 1:n_reps
    
    % call mfsobol.m with just the high-fidelity model to get Monte
    % Carlo estimate
    [sm,st,mu,sigsq] = mfsobol(fcns(1),d,w(1),stats,budget,vec(1),method);
    avg(n,1) = mu; 
    vr(n,1) = sigsq;
    mc_sm(n,:) = sm;
    mc_st(n,:) = st;
    
    % call mfsobol.m with full array of functions to get multifidelity
    % estimates
    [sm,st,mu,sigsq] = mfsobol(fcns,d,w,stats,budget,vec,method);
    avg(n,2) = mu; 
    vr(n,2) = sigsq;
    mf_sm(n,:) = sm;
    mf_st(n,:) = st;
end

%% PLOT ESTIMATOR SPREAD
warning('off','MATLAB:legend:IgnoringExtraEntries')
blue = [0       0.4470 0.7410];
red  = [0.8500  0.3250 0.0908];

% plot main effect sensitivity indices
figure(1); clf
h = boxplot([mc_sm(:,1), mf_sm(:,1), mc_sm(:,2), mf_sm(:,2),mc_sm(:,3),mf_sm(:,3),mc_sm(:,4),mf_sm(:,4),mc_sm(:,5),mf_sm(:,5),...
    mc_sm(:,6),mf_sm(:,6),mc_sm(:,7),mf_sm(:,7),mc_sm(:,8),mf_sm(:,8)],...
    'Colors',[blue; red; blue; red;blue;red;blue;red;blue;red;blue;red;blue;red;blue;red],'Whisker',10,...
    'labels',{'MC $s_m^1$','MF $s_m^1$','MC $s_m^2$','MF $s_m^2$','MC $s_m^3$','MF $s_m^3$',...
    'MC $s_m^4$','MF $s_m^4$','MC $s_m^5$','MF $s_m^5$','MC $s_m^6$','MF $s_m^6$',...
    'MC $s_m^7$','MF $s_m^7$','MC $s_m^8$','MF $s_m^8$'});
set(h,{'linew'},{2}); grid on
legend(flipud(findall(gca,'Tag','Box')), {'High-fidelity estimate','Multifidelity estimate'},...
    'Location','SouthWest','interpreter','latex'); legend boxoff
bp = gca;   bp.XAxis.TickLabelInterpreter = 'latex';
title([method,' main sensitivity estimates for Borehole model'],'interpreter','latex')

if strcmp(method,'Owen') || strcmp(method,'Saltelli')
    % plot total effect sensitivity indices
    figure(2); clf
    h = boxplot([mc_sm(:,1), mf_sm(:,1), mc_sm(:,2), mf_sm(:,2),mc_sm(:,3),mf_sm(:,3),mc_sm(:,4),mf_sm(:,4),mc_sm(:,5),mf_sm(:,5),...
    mc_sm(:,6),mf_sm(:,6),mc_sm(:,7),mf_sm(:,7),mc_sm(:,8),mf_sm(:,8)],...
        'Colors',[blue; red; blue; red;blue;red;blue;red;blue;red;blue;red;blue;red;blue;red],'Whisker',10,...
        'labels',{'MC $s_m^1$','MF $s_m^1$','MC $s_m^2$','MF $s_m^2$','MC $s_m^3$','MF $s_m^3$',...
    'MC $s_m^4$','MF $s_m^4$','MC $s_m^5$','MF $s_m^5$','MC $s_m^6$','MF $s_m^6$',...
    'MC $s_m^7$','MF $s_m^7$','MC $s_m^8$','MF $s_m^8$'});
    set(h,{'linew'},{2}); grid on
    hLegend = legend(flipud(findall(gca,'Tag','Box')), {'High-fidelity estimate','Multifidelity estimate'},...
        'Location','NorthEast','interpreter','latex'); legend boxoff
    bp = gca;   bp.XAxis.TickLabelInterpreter = 'latex';
    title([method,' total sensitivity estimates for Borehole model'],'interpreter','latex')
end

% plot variance estimates
figure(3); clf
histogram(vr(:,1),12,'facecolor',blue); hold on
histogram(vr(:,2),12,'facecolor',red,'facealpha',1)
legend({'High-fidelity','Multifidelity'},'interpreter','latex'); legend boxoff
title('Variance estimates for Borehole model','interpreter','latex')