function F2 = makeinterp(p)

Z1 = norminv(p, 0.10, 0.0161812);
Z2 = logninv(p, 7.71, 1.0056);
Z3 = (p * (115600 - 63070)) +63070 ;
Z4 = (p * (1110 - 990)) +990 ;
Z5 = (p * (116 - 63.1)) +63.1 ;
Z6 = (p * (820 - 700)) +700 ;
Z7 = (p * (1680 - 1120)) +1120 ;
Z8 = (p * (12045 - 9855)) +9855 ;

[a,b,c,d,e,f,g,h] = ndgrid(Z1, Z2, Z3, Z4, Z5, Z6, Z7, Z8);

term1 = 2*pi*c.*(d-f) ;
term2 = log(b./a).*(1 + ((2*g.*c)./(log(b./a).*(a.^2).*h)) + (c./e));

val = term1./term2 ;

F2 = griddedInterpolant(a,b,c,d,e,f,g,h, val, 'linear') ;

end