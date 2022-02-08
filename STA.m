function [T,width,omega,fidelity,psi_00,psi_11,psi,da,A,B] = STA(x,g,tf,wf,dt1,dt)
syms t a

dx = x(2)-x(1); M = (length(x))/2;n= 2*M;
dk=2*pi/(n*dx); 
k =(-n/2:n/2-1)*dk;%Define k-space grid

%%
w_i = 1;w_f = wf;
V_i = 0.5*w_i^2.*x.^2; N=1;

if g>100
    a_i = (1/w_i)^(2/3);% harmonic length
    
    psi_0 =  sqrt(1/a_i).*sech(x./a_i);
else
    a_i = sqrt(1/w_i);% harmonic length
    psi_0 = sqrt(N/a_i)*(1/pi)^(1/4)*exp(-x.^2./(2*a_i^2)); % Initial ansatz wavefunction
end

[psi_00,mu] = get_ground_state(psi_0,dt1,g,x,k,1,V_i,N); % actually wave-function
[E_k0,E_int0,E_V0,A,B,mu0,E_tot0] = get_energy(x,psi_00,g,V_i);


%% final states
V_f = 0.5*w_f^2.*x.^2; N=1;
if g>100
     a_f = (1/w_f)^(2/3);% harmonic length
    psi_1 =  sqrt(1/a_f).*sech(x./a_f);
else
    a_f = sqrt(1/w_f);% harmonic length
    psi_1= sqrt(N/a_f)*(1/pi)^(1/4)*exp(-x.^2./(2*a_f^2)); % Initial ansatz wavefunction
end
[psi_11,mu] = get_ground_state(psi_1,dt1,g,x,k,1,V_f,N); % actually wave-function
%% sta polynomial
ai=double(vpasolve(a*w_i^2-A/(a^3)-B/(a^2)==0,a,[0,inf]));%exact initial width
af=double(vpasolve(a*w_f^2-A/(a^3)-B/(a^2)==0,a,[0,inf]));%exact final width
a_i= max(ai); a_f = max(af);
T = 0:dt:tf;
a(t) = a_i- 6.*(a_i-a_f).*(t/tf).^5+15*(a_i-a_f).*(t/tf).^4 -10*(a_i-a_f).*(t/tf).^3;
da(t) = diff(a(t),1);
dda(t) = diff(a(t),2);
w2(t) = (-dda(t)+(A./a(t).^3+B./a(t).^2))./a(t);

%% time-propagation

Nt = tf/dt;
psi = psi_00;
ww2 = double(w2(T));
for i=1:Nt
    %t = dt*i;
    W2 = double(ww2(i));
    V = 0.5.*x.^2.*W2;
    [ psi ] = FFT( psi,V,g,dt,k,x);
    F = abs(sum(conj(psi).*psi_11).*dx).^2;
    
    %plot(x,abs(psi_00),x,abs(psi_11),'k:',x,abs(psi),'r--')
    %drawnow
end


fidelity = F;
omega = ww2;
width = double(a(T));
da = double(da(T));

end

