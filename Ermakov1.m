function dydt = Ermakov1(y,A,B,u)

dydt = zeros(2,1);
dydt(1) = y(2);%dot(a)
dydt(2) = A./y(1).^3+B./y(1).^2-u.*y(1);%


end

