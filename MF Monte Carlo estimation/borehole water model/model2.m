function f2 = model2(Z)
% first low-fidelity model for Ishigami function example
% f2 = sin(z1) + 0.95*a*sin(z2)^2 + b*z3^4*sin(z1)
%
% INPUTS
% Z         N-by-3 matrix of uncertain parameters, distributed ~U[-pi,pi]
% a,b       (optional) Ishigami function parameters, default a = 5, b = 0.1
%
% OUTPUT
% f2        N-by-1 vector of model evaluations at the uncertain inputs in Z
%
% AUTHOR
% Elizabeth Qian (elizqian@mit.edu) 14 June 2019




rw = Z(:,1);
r  = Z(:,2);
Tu = Z(:,3);
Hu = Z(:,4);
Tl = Z(:,5);
Hl = Z(:,6);
L  = Z(:,7);
Kw = Z(:,8);

frac1 = 2 * pi * Tu.*(Hu-Hl);

frac2a = 2*L.*Tu ./ ((rw.^2) .*Kw);
frac2b = Tu ./ Tl;
frac2 = log(r./rw) .* (1+frac2a+frac2b);


f2 = frac1 ./ frac2;