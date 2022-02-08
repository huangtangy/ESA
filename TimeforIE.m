function [ tf,w_min ] = TimeforIE(tf0,gN,k)
    syms a t
    w_0 = 250*2*pi;w_i = 250*2*pi/w_0;w_f = 2.5*2*pi/w_0;%intial and final frequency 
    a_i=double(vpasolve(a*w_i^2-1/(a^3)-gN/(a^2*sqrt(pi*2))==0,a,[0,inf]));
    a_f=double(vpasolve(a*w_f^2-1/(a^3)-gN/(a^2*sqrt(pi*2))==0,a,[0,inf]));
    beta = a_i;gamma =a_f;
    %psi_0=sqrt(N/a_i)*(1/pi)^(1/4)*exp(-x.^2/(2*a_i^2)); %Initial wavefunction
    %psi_T=sqrt(N/a_f)*(1/pi)^(1/4)*exp(-x.^2/(2*a_f^2)); %final wavefunction 
    Nt=100;w_min=1;j =1;
    
    while abs(w_min+1) > 0.01
        tf3 = tf0 - j*(0.02);
        at(t)= beta-6*(beta-gamma).*(t/tf3).^5+15*(beta-gamma).*(t/tf3).^4-10.*(beta-gamma).*(t/tf3).^3;%polynomial 
        %a = beta-6*(beta-gamma).*(T/tf3).^5+15*(beta-gamma).*(T/tf3).^4-10.*(beta-gamma).*(T/tf3).^3;%polynomial 
        a2(t) = diff(at,t,2); % ddot(a)
        T = 0:tf3/Nt:tf3;
        w2 = 1./at(T).^4+gN./(sqrt(2*pi).*at(T).^3)-a2(T)./at(T);%omega^2
        w_min = min(double(w2));
        %plot(T,w2)
        %drawnow
        if j>k
            'no solution found'
            tf = 0;
            w_min = 0;
            break
        end
        j = j+1;
        tf = tf3;
    end

    
end