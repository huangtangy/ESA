function [a_shell,Time,Aa_1,T,Tt1,Aa1,Ek1,Eint1,A0,B0] = quench_omega(wf,tf,g,dt,ii)
%mu_shell,Ev_shell,Eint_shell,Ek_shell,
%% calculating parameter A and B

syms t
m=1;%Dimensionless and defined the coordinate and kintic space
n=2^15;dx=0.01;x=(-n/2:n/2-1)*dx;
dk=2*pi/(n*dx); M =n/2;Nx = 2*M;
k =(-n/2:n/2-1)*dk;%Define k-space grid


Ek_shell = zeros();Ev_shell = zeros();Eint_shell=zeros();
a_shell = zeros();Time = zeros();mu_shell=zeros();
Etot_shell = zeros();
N =1;Nt = tf/dt;

g_i = g;omegaf =wf;
%% theory
omega_i = 1;a_i = sqrt(1/omega_i);a = a_i;
V = 0.5*m*omega_i.*x.^2;
psi_i=sqrt(N/a_i)*(1/pi)^(1/4)*exp(-x.^2/(2*a_i^2)); %Initial ansatz wavefunction
[psi_0,mu] = get_ground_state(psi_i,dt,g_i,x,k,m,V,N); % actually wave-function
[E_k,E_int,E_V,A0,B0,mu0,E_tot] = get_energy(x,psi_0,g,V);
[T,y] = ode45(@(t,y)Ermakov1(y,A0,B0,omegaf),[0:dt:tf],[a_i;0]);
aa= y(:,1);aa = aa./a_i;

%%
psi = psi_0;

for i = 1:1:Nt+1
    t = (i-1)*dt;
    
if t == 0
    V = 0.5*m*omega_i.*x.^2; 
    
else
    omega = omegaf;
    V = 0.5*m*omega^2.*x.^2; 
    %Omega = 0;
end

[ psi ] = FFT( psi,V,g,dt,k,x);

[E_k,E_int,E_V,A,B,mu,E_tot] = get_energy(x,psi,g,V);

mu_shell(i,1) = mu;
Ek_shell(i,1) = E_k;
Ev_shell(i,1) = E_V;
Eint_shell(i,1)=E_int;
Etot_shell(i,1) = E_k+E_V+E_int;
a = width(psi,x,dx);

a_shell(i,1) = a;
Time(i,1) =t;
end
a_shell = a_shell./a_shell(1);



Aa1=zeros();Tt1 = zeros();Aa_1=zeros();Tt_1 = zeros();Nt = length(Time)-1;
Ek1 = zeros();Ev1 = zeros();Eint1 = zeros();

for i = 1:1:ii+1
    if i == 1
        Aa1(i,1)= a_shell(1);
        Ek1(i,1)= Ek_shell(1);
        Ev1(i,1)= Ev_shell(1);
        Eint1(i,1)= Eint_shell(1);
        Aa_1(i,1)= aa(1);
        Tt1(i,1)= Time(1);
        
    else
        Aa1(i,1)= a_shell((i-1)*Nt/ii);
        Aa_1(i,1)= aa((i-1)*Nt/ii);
        Ek1(i,1)= Ek_shell((i-1)*Nt/ii);
        Ev1(i,1)= Ev_shell((i-1)*Nt/ii);
        Eint1(i,1)= Eint_shell((i-1)*Nt/ii);
        Tt1(i,1)= Time((i-1)*Nt/ii);
        
        %Aa1(i,1)= a_shell0((i-1)*Nt/10);
        %Tt1(i,1)= Time0((i-1)*Nt/10);
    end
end


end

