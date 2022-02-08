clc;close 
clf;clear
syms t
hbar=1;m=1;%Dimensionless and defined the coordinate and kintic space
n=2^11;dx=0.01;x=(-n/2:n/2-1)*dx;
dk=2*pi/(n*dx); M =n/2;Nx = 2*M;
k =(-n/2:n/2-1)*dk;%Define k-space grid
%dt = 0.0001;%time grid
%% ansatz of wave function
%% effective potential wf= 0
%{
g = 0.01; dt = 0.01;wf=0;tf=20;ii = 100;
[a_shell1,Time1,Aa_1,T1,Tt1,Aa1,Ek1,Eint1,A1,B1] = quench_omega(wf,tf,g,dt,ii);
g = 1; 
[a_shell2,Time2,Aa_2,T2,Tt2,Aa2,Ek2,Eint2,A2,B2] = quench_omega(wf,tf,g,dt,ii);
g = 100; dt= 0.001;
[a_shell3,Time3,Aa_3,T3,Tt3,Aa3,Ek3,Eint3,A3,B3] = quench_omega(wf,tf,g,dt,ii);

w_i =1;
%U1  = (0.5*(A1)./a_shell1.^2+B1./a_shell1)+0.5.*a_shell1.^2.*wf^2./w_i^2;
U_1 = (0.5*(A1)./Aa_1.^2+B1./Aa_1)+0.5.*Aa_1.^2.*wf^2./w_i^2;

%U2  = (0.5*(A2)./a_shell2.^2+B2./a_shell2)+0.5.*a_shell2.^2.*wf^2./w_i^2;
U_2 = (0.5*(A2)./Aa_2.^2+B2./Aa_2)+0.5.*Aa_2.^2.*wf^2./w_i^2;

%U3  = (0.5*(A3)./a_shell3.^2+B3./a_shell3)+0.5.*a_shell3.^2.*wf^2./w_i^2;
U_3 = (0.5*(A3)./Aa_3.^2+B3./Aa_3)+0.5.*Aa_3.^2.*wf^2./w_i^2;

%%
g = 0.01; dt = 0.01;wf=0.1;tf=100;ii = 100;
[a_shell11,Time11,Aa_11,T11,Tt11,Aa11,Ek11,Eint11,A11,B11] = quench_omega(wf,tf,g,dt,ii);
g = 1; 
[a_shell22,Time22,Aa_22,T22,Tt22,Aa22,Ek22,Eint22,A22,B22] = quench_omega(wf,tf,g,dt,ii);
g = 100; dt= 0.001;
[a_shell33,Time33,Aa_33,T33,Tt33,Aa33,Ek33,Eint33,A33,B33] = quench_omega(wf,tf,g,dt,ii);

w_i =1;
%U1  = (0.5*(A1)./a_shell1.^2+B1./a_shell1)+0.5.*a_shell1.^2.*wf^2./w_i^2;
U_11 = (0.5*(A11)./Aa_11.^2+B11./Aa_11)+0.5.*Aa_11.^2.*wf^2./w_i^2;

%U2  = (0.5*(A2)./a_shell2.^2+B2./a_shell2)+0.5.*a_shell2.^2.*wf^2./w_i^2;
U_22 = (0.5*(A22)./Aa_22.^2+B22./Aa_22)+0.5.*Aa_22.^2.*wf^2./w_i^2;

%U3  = (0.5*(A3)./a_shell3.^2+B3./a_shell3)+0.5.*a_shell3.^2.*wf^2./w_i^2;
U_33 = (0.5*(A33)./Aa_33.^2+B33./Aa_33)+0.5.*Aa_33.^2.*wf^2./w_i^2;


%%
subplot(1,2,2)
h1 = plot(a_shell1,U1./U1(1),'-r','Linewidth',1.5);
hold on
h11 = plot(a_shell2,U2./U2(1),'-k','Linewidth',1.5);
hold on 
h111 = plot(a_shell3,U3./U3(1),'-b','Linewidth',2.0);
hold on
legend([h1(1),h11(1),h111(1)],'g =0.01','g =1','g =100')
set(gca,'LineWidth',1.1,'FontSize',27,'Fontname','Times New Roman');
xlabel('$a$','interpret','latex')
ylabel('$\mathcal{U}(a)/\mathcal{U}(0)$','interpret','latex')

subplot(1,2,2)
H1 = plot(a_shell11,U11./U11(1),'-r','Linewidth',1.5);
hold on
H11 = plot(a_shell22,U22./U22(1),'-k','Linewidth',1.5);
hold on 
H111 = plot(a_shell33,U33./U33(1),'-b','Linewidth',2.0);
hold on

%{
plot(Aa_1(1:length(Aa_1)-1),diff(U_1),'r--','MarkerSize',10);
hold on
plot(Aa_2(1:length(Aa_2)-1),diff(U_2),'k-.','MarkerSize',10);
hold on
plot(Aa_3(1:length(Aa_3)-1),diff(U_3),'b:','MarkerSize',10);
hold on
%}
legend([H1(1),H11(1),H111(1)],'g =0.01','g =1','g =100')
set(gca,'LineWidth',1.1,'FontSize',27,'Fontname','Times New Roman');
xlabel('$a$','interpret','latex')
ylabel('$\mathcal{U}(a)/\mathcal{U}(0)$','interpret','latex')
%}

%%
shell_A = zeros();shell_B = zeros();shell_g = zeros();
UU= zeros();AA=zeros(); UUt = zeros();Time = zeros();
A_A =zeros();B_B =zeros();

%% w_f = 0.1;
wf = 0; 
for j = 1:1:6
G = [0.01,0.1,1,2,4,100];
g = G(j);dt = 0.001;tf = 20;

Nt = tf/dt;N =1;
w_i = 1;V_i = 0.5*m*w_i^2.*x.^2; 
a_i = sqrt(1/w_i);% harmonic length
psi_0=sqrt(N/a_i)*(1/pi)^(1/4)*exp(-x.^2./(2*a_i^2)); % Initial ansatz wavefunction
N0 = trapz(x,conj(psi_0).*psi_0);
[psi_0,mu] = get_ground_state(psi_0,dt,g,x,k,m,V_i,N); % actually wave-function
[E_k0,E_int0,E_V0,A,B,mu0,E_tot0] = get_energy(x,psi_0,g,V_i);

A_A(j,1) = A;B_B(j,1) = B;
end
%%

for kk=1:1:6
    %length(G)
    kk
A =A_A(kk);B= B_B(kk);
[T,y] = ode45(@(t,y)Ermakov1(y,A,B,wf),[0:dt:tf],[1;0]);
aa= y(:,1);da =(aa(end)-aa(1))/100;

a_T = aa;%1:200;%aa(1):da:aa(end);
U  = (0.5*(A)./a_T.^2+B./a_T)+0.5.*a_T.^2.*wf^2./w_i^2;

subplot(1,2,1)
plot(a_T,U./U(1),'Linewidth',3)
hold on
set(gca,'LineWidth',1.1,'FontSize',30,'Fontname','Times New Roman');
xlabel('$a$','interpret','latex')
ylabel('$\mathcal{U}(a)/\mathcal{U}(0)$','interpret','latex')
end
%% w_f = 0;
wf = 0.1; 
for j = 1:1:6
G = [0.01,0.1,1,2,4,100];
g = G(j);dt = 0.001;tf = 40;
Nt = tf/dt;N =1;
w_i = 1;V_i = 0.5*m*w_i^2.*x.^2; 
a_i = sqrt(1/w_i);% harmonic length
psi_0=sqrt(N/a_i)*(1/pi)^(1/4)*exp(-x.^2./(2*a_i^2)); % Initial ansatz wavefunction
N0 = trapz(x,conj(psi_0).*psi_0);
[psi_0,mu] = get_ground_state(psi_0,dt,g,x,k,m,V_i,N); % actually wave-function
[E_k0,E_int0,E_V0,A,B,mu0,E_tot0] = get_energy(x,psi_0,g,V_i);

A_A(j,1) = A;B_B(j,1) = B;
end
%%

for kk=1:1:6
    %length(G)
    kk
A =A_A(kk);B= B_B(kk);
[T,y] = ode45(@(t,y)Ermakov1(y,A,B,wf),[0:dt:tf],[1;0]);
aa= y(:,1);

da =(aa(end)-aa(1))/100;
a_T =aa;% 1:0.1:20;%aa(1):da:aa(end);
U  = (0.5*(A)./a_T.^2+B./a_T)+0.5.*a_T.^2.*wf^2./w_i^2;
subplot(1,2,2)
plot(a_T,U./U(1),'Linewidth',3)
hold on
set(gca,'LineWidth',1.1,'FontSize',30,'Fontname','Times New Roman');
xlabel('$a$','interpret','latex')
ylabel('$\mathcal{U}(a)/\mathcal{U}(0)$','interpret','latex')
end
