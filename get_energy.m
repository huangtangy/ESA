function [E_k,E_int,E_V,A,B,mu,E_tot] = get_energy(x,psi1,g,V)
%% get every part of instantaneous energy
dx = x(2)-x(1); M = (length(x))/2;
g_i =g;

n01 = (abs(psi1).^2); % density
%phi1 = asin(imag(psi1./sqrt(n01)));
%phi2 = asin(imag(psi2./sqrt(n02)));
x1=((-M+1):1:M-1).*dx; 

E_k = 0.5*trapz(x1,(diff(sqrt(n01),1)./dx).^2);
E_int = g_i*trapz(x,n01.^2);
E_V = trapz(x,V.*n01);


dd1 = diff(sqrt(n01),2)./dx^2;
G_0 = g_i*(n01(M));
sn01 = sqrt(n01);
D_0 = -0.5*(dd1(M))/(sn01(M));

mu = G_0+D_0;% chemical potential
E_tot = E_int+E_k+E_V;

B = (G_0-E_int)/E_V;
A = (D_0-E_k)/E_V;


end

