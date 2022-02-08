function [a] = width(psi1,x,dx)
a = sqrt(2*(trapz(x,conj(psi1).*x.^2.*psi1)-(trapz(x,conj(psi1).*x.*psi1))^2));
a = abs(a);
end

