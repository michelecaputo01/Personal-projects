function f1 = model1(Z)

a = 5 ; 
b = 0.1 ;

f1 = sin(Z(:,1)) +  a*sin(Z(:,2)).^2 + b*Z(:,3).^4.*sin(Z(:,1));