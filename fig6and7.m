clc;close 
clf;clear
syms t
hbar=1;m=1;%Dimensionless and defined the coordinate and kintic space
n=2^14;dx=0.01;
%n=2^13;dx=0.1;
x=(-n/2:n/2-1)*dx;
dk=2*pi/(n*dx); M =n/2;Nx = 2*M;
k =(-n/2:n/2-1)*dk;%Define k-space grid
%dt = 0.0001;%time grid
%%

TT = zeros();AA = zeros();WW2= zeros();Psi=zeros();Psi1=zeros();
EE = zeros()
for i = 1:1:3
    i
G = [0.01,1,100];
tf = 1;wf = 0.1;dt1 = 0.001;dt = 0.01;
g = G(i);
[T,width,omega,fidelity,psi_00,psi_11,psi,da,A,B] = STA(x,g,tf,wf,dt1,dt);
a = width; w2 = omega;% omega^2

for j = 1:length(T)
    TT(i,j) = T(j);
    AA(i,j) = a(j);
    WW2(i,j) = w2(j);
    EE(i,j) = 0.5*da(j)+A/(2*a(j)^2)+B/a(j)+0.5*(w2(j))*a(j)^2;
end

for j = 1:length(psi)
    Psi(i,j) = psi(j);
    Psi1(i,j) = psi_11(j);
end

end



figure(1)
subplot(2,2,1)
h1 = plot(TT(1,:),WW2(1,:),'r--',TT(2,:),WW2(2,:),'b-',TT(3,:),WW2(3,:),'k:','Linewidth',1);
H1 = legend([h1(1),h1(2),h1(3)],'$$g = 0.01$$','$$g = 1$$','$$g = 100$$');
set(H1,'Interpreter','latex','Location','SouthOutside') %设置legend为latex解释器显示分式
set(gca,'LineWidth',1.1,'FontSize',27,'Fontname','Times New Roman');
xlabel('$t/t_f$','interpret','latex')
ylabel('$\omega^2(t)$','interpret','latex')

subplot(2,2,2)
h2 = plot(TT(1,:),AA(1,:),'r--',TT(2,:),AA(2,:),'b-',TT(3,:),AA(3,:),'k:','Linewidth',1);
H2 = legend([h2(1),h2(2),h2(3)],'$$g = 0.01$$','$$g = 1$$','$$g = 100$$');
set(H2,'Interpreter','latex','Location','SouthOutside') %设置legend为latex解释器显示分式
set(gca,'LineWidth',1.1,'FontSize',27,'Fontname','Times New Roman');
xlabel('$t/t_f$','interpret','latex')
ylabel('$a(t)/a(0)$','interpret','latex')

subplot(2,1,2)
h33 = plot(x,abs(Psi(1,:)),'c-',x,abs(Psi(2,:)),'c-',x,abs(Psi(3,:)),'c-','Linewidth',2);
hold on
h3 = plot(x,abs(Psi1(1,:)),'r--',x,abs(Psi1(2,:)),'b-.',x,abs(Psi1(3,:)),'k:','Linewidth',2);
hold on
H3 = legend([h3(1),h3(2),h3(3)],'$$g = 0.01$$','$$g = 1$$','$$g = 100$$');
set(H3,'Interpreter','latex','Location','SouthOutside') %设置legend为latex解释器显示分式
set(gca,'LineWidth',1.1,'FontSize',27,'Fontname','Times New Roman');
xlabel('$x$','interpret','latex')
ylabel('$\tilde{\psi}_f(x)/\psi(x,t_f)$','interpret','latex')

%}
%%


load('FF1.mat')
load('GG1.mat')

figure(2)
plot(GG,FF2,GG,FF1,GG,FF3,'Linewidth',1)
hold on
h = legend('$$\omega = 0.2$$','$$\omega = 0.1$$','$$\omega = 0.05$$');
set(h,'Interpreter','latex','Location','SouthOutside') %设置legend为latex解释器显示分式
set(gca,'LineWidth',1.1,'FontSize',30,'Fontname','Times New Roman');
xlabel('$g$','interpret','latex')
ylabel('$F$','interpret','latex')




figure(3)
%h2 = plot(AA(1,:),EE(1,:),'r--',AA(2,:),EE(2,:),'b-',AA(3,:),EE(3,:),'k:','Linewidth',1);
h2 = plot(TT(1,:),EE(1,:),'r--',TT(2,:),EE(2,:),'b-',TT(3,:),EE(3,:),'k:','Linewidth',1);
H2 = legend([h2(1),h2(2),h2(3)],'$$g = 0.01$$','$$g = 1$$','$$g = 100$$');
set(H2,'Interpreter','latex','Location','SouthOutside') %设置legend为latex解释器显示分式
set(gca,'LineWidth',1.1,'FontSize',24,'Fontname','Times New Roman');
xlabel('$t/t_f$','interpret','latex')
%xlabel('$a(t)/a(0)$','interpret','latex')
ylabel('$E(a)$','interpret','latex')

