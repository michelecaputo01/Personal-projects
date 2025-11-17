function [a, P, x_MLM] = MLM_unmixing(x, E, P_min)

% Perform unmixing according to MLM; x is the spectral matrix (with shape n_bands 
% x lines*columns); E the endmember matrix (with shape n_bands x n_endmembers);
% a is the abundance matrix (with shape n_endmembers x lines*columns); P is
% the vector containing the P values (of the MLM) for each pixel (with length
% lines*columns); P_min is the minimal allowed value for P

% Mostly extracted from 'Demo on the MLM model'. Here follows a list of the
% modifications of the original code:

% - The original objective function has been divided by the number of
%   included bands, to avoid penalizing the hyperspectral images with many
%   bands
% - The script has been converted into a MATLAB function
% - Only LMM and MLM have been included
% - A counter of the iterations has been added
% - Some comments have been removed

% Reference:
% R. Heylen, P. Scheunders, "A multilinear mixing model for nonlinear 
% spectral unmixing," IEEE Tran. Geosci. Remote Sens., vol. PP, no. 99., 
% pp. 1-12, 2015
%
% Rob Heylen, University of Antwerp, 2015



%% Optimizer settings

options = optimset('fmincon');
options = optimset(options,'Display','off','Algorithm','sqp','MaxFunEvals',50000,...
                   'TolFun',1e-10,'TolCon',1e-8,'GradObj','off');

[d,N]=size(x);
[~,p]=size(E);



%% Linear unmixing to initialize the MLM solution

a_LMM=zeros(p,N);
x_LMM=zeros(d,N);

for i=1:N
    
    if(mod(i,500)==0)
        disp(i);
    end
    
    a_ini=ones(p,1)/p; % Initialize with mean of endmembers

    a_LMM(:,i) = fmincon(@(a) sum((x(:,i)-model_LMM(a,E)).^2)/d, a_ini,[],...
                         [],ones(1,p),1,zeros(p,1),ones(p,1),[],options);
    x_LMM(:,i) = model_LMM(a_LMM(:,i),E);
end



%% MLM unmixing

a_MLM=zeros(p+1,N); % The first p variables are the abundances, the p+1'th variable contains the P value
x_MLM=zeros(d,N);

for i=1:N
    
    if(mod(i,500)==0)
        disp(i);
    end
    
    a_ini=[a_LMM(:,i); 0.0]; % Initialize with linear unmixing results, P=0
                             % Sum-to-one applies to abundances, not P. P is restricted to [P_min,1]
                             
    % I replaced the SE with the MSE, in the objective function
    a_MLM(:,i) = fmincon(@(a) sum((x(:,i)-model_MLM(a,E)).^2)/d, a_ini,[],...
                         [],[ones(1,p) 0],1,[zeros(p,1); P_min],ones(p+1,1),[],options);
    x_MLM(:,i) = model_MLM(a_MLM(:,i),E);
end



%% Returned matrices

a = a_MLM(1:p,:);
P = a_MLM(p+1,:);