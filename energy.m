function [EE] = energy(a,A,B,w2)
EE = zeros();
da = diff(a,1);
for j = 1:length(a)-1
    EE(j) = 0.5*da(j)+A/(2*a(j)^2)+B/a(j)+0.5*w2*a(j)^2;
end
end

