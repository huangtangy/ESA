function [LA] = period_driving(x,k,g,tf,Wt)
%% initial states
dx = x(2)-x(1);
w_i = 1;
V_i = 0.5*w_i^2.*x.^2; N=1;

if g>100
    a_i = (1/w_i)^(2/3);% harmonic length
    
    psi_0 =  sqrt(1/a_i).*sech(x./a_i);
else
    a_i = sqrt(1/w_i);% harmonic length
    psi_0 = sqrt(N/a_i)*(1/pi)^(1/4)*exp(-x.^2./(2*a_i^2)); % Initial ansatz wavefunction
end
dt1 = 0.001;
[psi_00,mu] = get_ground_state(psi_0,dt1,g,x,k,1,V_i,N); % actually wave-function
[E_k0,E_int0,E_V0,A,B,mu0,E_tot0] = get_energy(x,psi_00,g,V_i);
A+B
%%
dt =0.1;LA=zeros();
Nt = tf/dt;
psi = psi_00;
for i =1:Nt
    
    %t = dt*i;
    W2 = double(Wt(i));
    V = 0.5.*x.^2.*W2;
    [ psi ] = FFT( psi,V,g,dt,k,x);
    la = (sum(conj(psi).*psi_00).*dx); 
    LA(i,1) = la;
    %plot(x,abs(psi_00),x,abs(psi_11),'k:',x,abs(psi),'r--')
    %drawnow
end




end

