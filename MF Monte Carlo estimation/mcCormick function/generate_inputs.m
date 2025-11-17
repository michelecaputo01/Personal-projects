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

Z(:, 1) = (rand(N, 1) * (4 + 1.5)) - 1.5;
Z(:, 2) = (rand(N, 1) * (4 + 3)) -3 ;