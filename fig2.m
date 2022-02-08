clc;close 
clf;clear
syms t
hbar=1;m=1;%Dimensionless and defined the coordinate and kintic space
n=2^15;dx=0.1;x=(-n/2:n/2-1)*dx;
dk=2*pi/(n*dx); M =n/2;Nx = 2*M;
k =(-n/2:n/2-1)*dk;%Define k-space grid
%dt = 0.0001;%time grid
%% g = inf
Dx1= zeros();T_T1 = zeros();Dx2= zeros();T_T2 = zeros();
WW = zeros();Dx22= zeros();
tf = 50;dt = 0.001;wi = 1;
for i = 1:1:99
i
wf = 0.01+(i-1)*0.01;
if wf< 0.06
    tf = 200;
else
    tf = 50;
end
%% g = 0
[T,y] = ode45(@(t,y)Ermakov1(y,1,0,wf),[0:dt:tf],[1;0]);
aa= y(:,1);
P = pi/wf;amp = (wi^2-wf^2)/wf^2;
aa1 = sqrt(1+(amp).*sin(pi.*T./P).^2);
a_min = min(aa); a_max = max(aa);

amp = sqrt(1+(wi^2-wf^2)/wf^2);
Dx1(i,1)= a_max-a_min;

TT = find_period1(aa,T);
T_T1(i,1) = TT;
%% g = inf
[T2,y] = ode45(@(t,y)Ermakov1(y,0,1,wf),[0:dt:tf],[1;0]);
aa_2= y(:,1);
%P2 = 0.89*P;amp2 = (1.337*amp);%sqrt((wi/wf)^3-1)%
%%

%P2 = exp(0.965*log(P)+0.1958);
P2 = exp((0.029*wf^3-0.07*wf^2+0.17*wf+1)*log(P));
amp2 = 2*(wi^2-wf^2)/(wf^2); %pridict
aa2 = sqrt(2+amp2.*sin(pi.*T2./P2).^2)-(sqrt(2)-1);

a_min = min(aa_2); a_max = max(aa_2);
a_min1 = min(aa2); a_max1 = max(aa2);

Dx2(i,1)= a_max-a_min;
Dx22(i,1)= a_max1-a_min1;
%%

TT1 = find_period1(aa_2,T2);
TTT1 = find_period1(aa2,T2);

T_T2(i,1) = TT1;
WW(i,1) = wf;

end
%{
figure(1)
plot(T,aa,T,aa1)
figure(2)
plot(T2,aa_2,T2,aa2)
figure(3)
plot(T,sin(T),T,sqrt(sin(T).^2),'--r')
%}

figure(1)
subplot(1,2,1)
plot(WW,log(T_T2)./log(T_T1),'b-','Linewidth',1.5)
hold on
x = WW;
plot(x,0.036.*x.^3 - 0.082.*x.^2 + 0.17.*x + 1,'r--','Linewidth',3)
h = legend('Numerical data','$$y = 0.036*x^{3} - 0.082*x^{2} + 0.17*x + 1$$');
set(h,'Interpreter','latex','Location','SouthOutside') %设置legend为latex解释器显示分式
set(gca,'LineWidth',1.1,'FontSize',30,'Fontname','Times New Roman');
xlabel('$\omega_f$','interpret','latex')
ylabel('$log(T_2)/log(T_1)$','interpret','latex')

subplot(1,2,2)
plot(Dx1,Dx2,'b-','Linewidth',1.0)
hold on
plot(Dx1,sqrt(2).*Dx1,'r--','Linewidth',1.5)
hold on
h = legend('Numerical data','$$y = \sqrt{2}x$$');
set(h,'Interpreter','latex','Location','SouthOutside') %设置legend为latex解释器显示分式
set(gca,'LineWidth',1.1,'FontSize',30,'Fontname','Times New Roman');
xlabel('$Ap_1$','interpret','latex')
ylabel('$Ap_2$','interpret','latex')


%}
%{
subplot(2,2,3)
plot(log(T_T1),log(T_T2),'b-','Linewidth',1.0)
hold on
plot(log(T_T1),0.965*log(T_T1)+0.1958,'r--','Linewidth',1.5)
set(gca,'LineWidth',1.1,'FontSize',30,'Fontname','Times New Roman');
xlabel('$log(T_1)$','interpret','latex')
ylabel('$log(T_2)$','interpret','latex')

subplot(2,2,4)
plot(WW,Dx2./Dx1,'b-','Linewidth',1.5)
hold on
set(gca,'LineWidth',1.1,'FontSize',30,'Fontname','Times New Roman');
xlabel('$\omega_f$','interpret','latex')
ylabel('$Ap_2/Ap_1$','interpret','latex')

%}