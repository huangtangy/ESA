clc;close 
clf;clear
syms t
hbar=1;m=1;%Dimensionless and defined the coordinate and kintic space
n=2^14;dx=0.01;x=(-n/2:n/2-1)*dx;
dk=2*pi/(n*dx); M =n/2;Nx = 2*M;k =(-n/2:n/2-1)*dk;%Define k-space grid
%dt = 0.0001;%time grid
%% free expansion
g = 0.01; dt = 0.01;wf=0; tf=20; ii = 20;
[a_shell11,Time11,Aa_11,T1,Tt11,Aa1,Ek1,Eint1,A11,B11] = quench_omega(wf,tf,g,dt,ii);
g = 1; 
[a_shell22,Time22,Aa_22,T2,Tt22,Aa2,Ek2,Eint2,A22,B22] = quench_omega(wf,tf,g,dt,ii);
g = 100; dt= 0.001;
[a_shell33,Time33,Aa_33,T3,Tt33,Aa3,Ek3,Eint3,A33,B33] = quench_omega(wf,tf,g,dt,ii);

subplot(2,1,1)
plot(Time11,a_shell11,'--k','Linewidth',1.0)
hold on
plot(Time22,a_shell22,'-.k','Linewidth',1.0)
hold on 
plot(Time33,a_shell33,':k','Linewidth',1.5)
hold on
h1 = plot(Tt11,Aa_11,'ro','MarkerSize',10);
hold on
h11 = plot(Tt22,Aa_22,'kd','MarkerSize',10);
hold on
h111 = plot(Tt33,Aa_33,'bs','MarkerSize',10);
hold on
legend([h1(1),h11(1),h111(1)],'g =0.01','g =1','g =100')
set(gca,'LineWidth',1.1,'FontSize',27,'Fontname','Times New Roman');
xlabel('$t$','interpret','latex')
ylabel('$a(t)/a(0)$','interpret','latex')

%% quench breath model wf = 0.1

g = 0.01; dt = 0.01;wf=0.1;tf=100;ii = 100;
[a_shell1,Time1,Aa_1,T1,Tt1,Aa1,Ek1,Eint1,A1,B1] = quench_omega(wf,tf,g,dt,ii);
g = 1; 
[a_shell2,Time2,Aa_2,T2,Tt2,Aa2,Ek2,Eint2,A2,B2] = quench_omega(wf,tf,g,dt,ii);
g = 100; dt= 0.001;
[a_shell3,Time3,Aa_3,T3,Tt3,Aa3,Ek3,Eint3,A3,B3] = quench_omega(wf,tf,g,dt,ii);

subplot(2,1,2)
plot(Time1,a_shell1,'--k','Linewidth',1.0)
hold on
plot(Time2,a_shell2,'-.k','Linewidth',1.0)
hold on 
plot(Time3,a_shell3,':k','Linewidth',1.5)
hold on
h1 = plot(Tt1,Aa_1,'ro','MarkerSize',10);
hold on
h11 = plot(Tt2,Aa_2,'kd','MarkerSize',10);
hold on
h111 = plot(Tt3,Aa_3,'bs','MarkerSize',10);
hold on
legend([h1(1),h11(1),h111(1)],'g =0.01','g =1','g =100')
set(gca,'LineWidth',1.1,'FontSize',27,'Fontname','Times New Roman');
xlabel('$t$','interpret','latex')
ylabel('$a(t)/a(0)$','interpret','latex')


%% energy

figure(2)

subplot(2,1,1)
h1 =plot(a_shell11(1:2000),energy(a_shell11,A11,B11,0.1),'--k','Linewidth',1.0);
hold on
h11 =plot(a_shell22(1:2000),energy(a_shell22,A22,B22,0.1),'-.k','Linewidth',1.0);
hold on 
h111 = plot(a_shell33(1:20000),energy(a_shell33,A33,B33,0.1),':k','Linewidth',1.5);
hold on
legend([h1(1),h11(1),h111(1)],'g =0.01','g =1','g =100')
set(gca,'LineWidth',1.1,'FontSize',27,'Fontname','Times New Roman');
ylabel('$E(a)$','interpret','latex')
xlabel('$a(t)/a(0)$','interpret','latex')

subplot(2,1,2)
%h1 =plot(a_shell1(1:10000),energy(a_shell1,A11,B11,0),'--k','Linewidth',0.5);
h1 =plot(Time1(1:10000),energy(a_shell1,A11,B11,0),'--k','Linewidth',1);
hold on
%h11 =plot(a_shell2(1:10000),energy(a_shell2,A22,B22,0),'-k','Linewidth',0.5);
h11 =plot(Time2(1:10000),energy(a_shell2,A22,B22,0),'-k','Linewidth',1);
hold on 
h111 = plot(Time3(1:100000),energy(a_shell3,A33,B33,0),':k','Linewidth',1.5);
%h111 = plot(a_shell3(1:100000),energy(a_shell3,A33,B33,0),':k','Linewidth',0.5);
hold on
legend([h1(1),h11(1),h111(1)],'g =0.01','g =1','g =100')
set(gca,'LineWidth',1.1,'FontSize',27,'Fontname','Times New Roman');
ylabel('$E(a)$','interpret','latex')
xlabel('$t/t_f$','interpret','latex')




%% effective potential
%{
w_i =1;
U1  = (0.5*(A1)./a_shell1.^2+B1./a_shell1)+0.5.*a_shell1.^2.*wf^2./w_i^2;
U_1 = (0.5*(A1)./Aa_1.^2+B1./Aa_1)+0.5.*Aa_1.^2.*wf^2./w_i^2;

U2  = (0.5*(A2)./a_shell2.^2+B2./a_shell2)+0.5.*a_shell2.^2.*wf^2./w_i^2;
U_2 = (0.5*(A2)./Aa_2.^2+B2./Aa_2)+0.5.*Aa_2.^2.*wf^2./w_i^2;

U3  = (0.5*(A3)./a_shell3.^2+B3./a_shell3)+0.5.*a_shell3.^2.*wf^2./w_i^2;
U_3 = (0.5*(A3)./Aa_3.^2+B3./Aa_3)+0.5.*Aa_3.^2.*wf^2./w_i^2;

%%

subplot(2,1,2)
h1 = plot(a_shell1,U1./U1(1),'-r','Linewidth',1.5);
hold on
h11 = plot(a_shell2,U2./U2(1),'-k','Linewidth',1.5);
hold on 
h111 = plot(a_shell3,U3./U3(1),'-b','Linewidth',2.0);
hold on
%{
plot(Aa_1(1:length(Aa_1)-1),diff(U_1),'r--','MarkerSize',10);
hold on
plot(Aa_2(1:length(Aa_2)-1),diff(U_2),'k-.','MarkerSize',10);
hold on
plot(Aa_3(1:length(Aa_3)-1),diff(U_3),'b:','MarkerSize',10);
hold on
%}
legend([h1(1),h11(1),h111(1)],'g =0.01','g =1','g =100')
set(gca,'LineWidth',1.1,'FontSize',27,'Fontname','Times New Roman');
xlabel('$a$','interpret','latex')
ylabel('$\mathcal{U}(a)/\mathcal{U}(0)$','interpret','latex')
%}