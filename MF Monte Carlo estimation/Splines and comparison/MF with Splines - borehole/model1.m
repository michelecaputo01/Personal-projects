function f1 = model1(Z)

rw = Z(:,1);
r = Z(:,2);
tu = Z(:,3);
hu = Z(:,4);
tl = Z(:,5);
hl = Z(:,6);
l = Z(:,7);
kw = Z(:,8);

term1 = 2*pi*tu.*(hu-hl) ;
term2 = log(r./rw).*(1 + ((2*l.*tu)./(log(r./rw).*(rw.^2).*kw)) + (tu./tl));

f1 = term1./term2 ;