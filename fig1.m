clc;close 
clf;clear
syms t
hbar=1;m=1;%Dimensionless and defined the coordinate and kintic space
n=2^16;dx=0.001;x=(-n/2:n/2-1)*dx;
dk=2*pi/(n*dx); M =n/2;Nx = 2*M;
k =(-n/2:n/2-1)*dk;%Define k-space grid
tic
%dt = 0.0001;%time grid
%% ansatz of wave function
%{
N = 1; %atom number
w_i = 1;V_i = 0.5*m*w_i^2.*x.^2; 
g_i = 10;%initial parameters frequency and interaction
a_i = sqrt(1/w_i);% harmonic length

psi_0=sqrt(N/a_i)*(1/pi)^(1/4)*exp(-x.^2/(2*a_i^2)); % Initial ansatz wavefunction
plot(x,abs(psi_0).^2,'b-')
hold on
[psi_0,mu] = get_ground_state(psi_0,dt,g_i,x,k,m,V_i,N); % actually wave-function
plot(x,abs(psi_0).^2,'r--')
%}
%% calculating parameter A and B


%%
shell_A1 = zeros();shell_B1 = zeros();shell_g1 = zeros();
dt = 0.0001;%time grid
Ni = 1;
for i = 1:1:Ni
    i
    
j = 3+(i-1)*(6/Ni);
N = 1;g =10^(j);
w_i = 1;a_i = sqrt(1/w_i);
a = a_i; 

g_i = g*N;
%x = x./a;dx= dx./a;%re-scale
%dk=2*pi/(n*dx);k =(-n/2:n/2-1)*dk;
w0= w_i;V_i = 0.5*m*w_i^2.*x.^2; %initial parameters frequency and interaction

psi_i=sqrt(N/a_i)*(1/pi)^(1/4)*exp(-x.^2/(2*a_i^2)); %Initial ansatz wavefunction
[psi_0,mu] = get_ground_state(psi_i,dt,g_i,x,k,m,V_i,N); % actually wave-function


n0 = (abs(psi_0).^2); % density
sig_s = trapz(x,x.^2.*n0); % sigma squre
E_int = 0.5*g_i*trapz(x,n0.^2); 
x1=((-M+1):1:M-1).*dx; 
E_k = 0.5*trapz(x1,(diff(sqrt(n0))./dx).^2);
E_V =trapz(x,V_i.*n0);
 
dd = diff(n0,2)./dx^2;
D_0 = -(0.25./(a^2.*n0(M))).*dd(M);

A = ((2/(sig_s*w0^2))*(D_0-E_k));
B = ((2/(sig_s*w0^2))*(g_i*n0(M)-2*E_int));

muu = g_i*n0(M)+D_0;%mu
AB = ((2/(sig_s*w0^2))*(muu-E_k-2*E_int));
%delta = abs(muu-E);

shell_A1(1,i) = A; 
shell_B1(1,i)= B;
shell_g1(1,i) = g_i;
A+B

end

%{
dt = 0.0001;%time grid
parfor i = 51:1:100
    i
    
j = -3+(i-1)*0.06;
N = 1;g =10^(j);
w_i =1;a_i = sqrt(1/w_i);
a = a_i;

g_i = g*N;
w0= w_i;V_i = 0.5*m*w_i^2.*x.^2; %initial parameters frequency and interaction


psi_i=sqrt(N/a_i)*(1/pi)^(1/4)*exp(-x.^2/(2*a_i^2)); %Initial ansatz wavefunction
[psi_0,mu] = get_ground_state(psi_i,dt,g_i,x,k,m,V_i,N); % actually wave-function

n0 = (abs(psi_0).^2); % density
sig_s = trapz(x,x.^2.*n0); % sigma squre
E_int = 0.5*g_i*trapz(x,n0.^2); 
x1=((-M+1):1:M-1).*dx; 
E_k = 0.5*trapz(x1,(diff(sqrt(n0))./dx).^2);
E_V =trapz(x,V_i.*n0);
 
dd = diff(n0,2)./dx^2;
D_0 = -(0.25./(a^2.*n0(M))).*dd(M);
A = ((2/(sig_s*w0^2))*(D_0-E_k));
B = ((2/(sig_s*w0^2))*(g_i*n0(M)-2*E_int));
muu = g_i*n0(M)+D_0;%mu
     
AB = ((2/(sig_s*w0^2))*(muu-E_k-2*E_int));
%delta = abs(muu-E);
A+B

shell_A1(1,i) = A; shell_B1(1,i)= B;
shell_g1(1,i) = g_i;

end


%}

subplot(2,1,1)
plot(shell_g1,shell_A1,'--r',shell_g1,shell_B1,'-b','Linewidth',1.0)
hold on
plot(shell_g1,shell_A1+shell_B1,'-.k','Linewidth',1.0)
set(gca,'LineWidth',1.1,'FontSize',20,'Fontname','Times New Roman');
xlabel('g')
ylabel('A/B')
toc