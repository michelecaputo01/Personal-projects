function F2 = makeinterp(p)

Z1 = (p * 2 * pi) - pi;
Z2 = (p * 2 * pi) - pi;
Z3 = (p * 2 * pi) - pi;

a = 5 ;
b = 0.1 ;

[x,y,z] = ndgrid(Z1, Z2, Z3);

val = sin(x) +  a*sin(y).^2 + b*(z.^4).*sin(x) ;

F2 = griddedInterpolant(x,y,z, val, 'linear') ;

end