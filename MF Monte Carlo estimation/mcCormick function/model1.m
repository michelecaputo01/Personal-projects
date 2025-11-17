function f1 = model1(Z)
% high-fidelity model for Ishigami function example
% f1 = sin(z1) + a*sin(z2)^2 + b*z3^4*sin(z1)
%
% INPUTS
% Z         N-by-2 matrix of uncertain parameters, Z1 U[-1.5,4], Z2 U[-3,4]
% a,b       (optional) Ishigami function parameters, default a = 5, b = 0.1
%
% OUTPUT
% f1        N-by-1 vector of model evaluations at the uncertain inputs in Z
%
% AUTHOR
% Elizabeth Qian (elizqian@mit.edu) 14 June 2019


f1 = sin(Z(:,1) + Z(:,2) ) + (Z(:,1) + Z(:,2) ).^2 -1.5.*Z(:,1) + 2.5.*Z(:,2) + 1 ;
