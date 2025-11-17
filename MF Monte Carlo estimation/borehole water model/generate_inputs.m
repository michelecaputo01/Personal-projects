function Z = generate_inputs(N)
% generates random inputs for mcCormick function example
%
% INPUT
% N     number of inputs to generate
%
% OUTPUT
% Z     N-by-2 matrix of random McCormick function inputs
%
% AUTHOR
% Elizabeth Qian (elizqian@mit.edu) 14 June 2019

Z(:,1) = normrnd(0.10,0.0161812,[N,1]); %radius of borehole(m) rw
Z(:,2) =  lognrnd(7.71, 1.0056, [N,1]); %radius of influence(m) r
Z(:,3)  = (rand(N, 1) * (115600 - 63070)) +63070 ; %transmissivity of upper aquifer (m2/yr) Tu
Z(:,4)  = (rand(N, 1) * (1110 - 990)) +990 ;  %potentiometric head of upper aquifer (m) Hu
Z(:,5)  = (rand(N, 1) * (116 - 63.1)) +63.1  ; %transmissivity of lower aquifer (m2/yr) Tl
Z(:,6)  = (rand(N, 1) * (820 - 700)) +700 ; %potentiometric head of lower aquifer (m) Hl
Z(:,7)  = (rand(N, 1) * (1680 - 1120)) +1120 ; %length of borehole(mm) L
Z(:,8)  = (rand(N, 1) * (12045 - 9855)) +9855; %hydraulic conductivity of borehole (m/yr) Kw
