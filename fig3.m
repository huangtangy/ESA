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
WW = zeros();Dx22= zeros();Dx1= zeros();Dx11= zeros();
tf = 50;dt = 0.001;wi = 1;Width = zeros();Width1 = zeros();Width2=zeros();
width = zeros();width1 = zeros();width2=zeros();

for i = 1:1:5
i
wf = 0.01;
if wf< 0.06
    tf = 500;
else
    tf = 50;
end

%% g= 1

G = [0.01,0.5,1,2,100];
g = G(i)
w_i = 1;V_i = 0.5*w_i^2.*x.^2; N=1;
a_i = sqrt(1/w_i);% harmonic length
psi_0=sqrt(N/a_i)*(1/pi)^(1/4)*exp(-x.^2./(2*a_i^2)); % Initial ansatz wavefunction

if g > 90
    dt1 = 0.0001; 
else
    dt1 = 0.001;
end

[psi_0,mu] = get_ground_state(psi_0,dt1,g,x,k,1,V_i,N); % actually wave-function
[E_k0,E_int0,E_V0,A,B,mu0,E_tot0] = get_energy(x,psi_0,g,V_i);

%% g = 0
[T,y] = ode45(@(t,y)Ermakov1(y,1,0,wf),[0:dt:tf],[1;0]);
aa= y(:,1);
P = pi/wf;amp = 1*(wi^2-wf^2)/wf^2;
aa_ = sqrt(1+(amp).*sin(pi.*T./P).^2);

%% g = inf
[T1,y1] = ode45(@(t,y)Ermakov1(y,0,1,wf),[0:dt:tf],[1;0]);
aa_1= y1(:,1);

%P1 = exp((0.029*wf^3-0.07*wf^2+0.17*wf+1)*log(P));
P1 = exp((- 0.031*wf^4+0.098*wf^3 - 0.12*wf^2 +0.18*wf+1)*log(P));
 
amp1 =1*2*(wi^2-wf^2)/(wf^2); %pridict
aa1 = sqrt(2+amp1.*sin(pi.*T1./P1).^2)-(sqrt(2)-1);

%%
[T2,y2] = ode45(@(t,y)Ermakov1(y,A,B,wf),[0:dt:tf],[1;0]);
aa_2= y2(:,1);

aa2 = (A*aa+B*aa1);
%%
figure(1)
subplot(3,1,1)
plot(T2,aa_2,'k-','Linewidth',1.0)
hold on
plot(T,aa2,'k:','Linewidth',1.5)
hold on
set(gca,'LineWidth',1.1,'FontSize',30,'Fontname','Times New Roman');
xlabel('$t$','interpret','latex')
ylabel('$a(t)$','interpret','latex')


end


for i = 1:1:5
i
wf = 0.1;
if wf< 0.06
    tf = 200;
else
    tf = 50;
end

%% g= 1

G = [0.01,0.5,1,2,100];
g = G(i)
w_i = 1;V_i = 0.5*w_i^2.*x.^2; N=1;
a_i = sqrt(1/w_i);% harmonic length
psi_0=sqrt(N/a_i)*(1/pi)^(1/4)*exp(-x.^2./(2*a_i^2)); % Initial ansatz wavefunction

if g > 90
    dt1 = 0.0001; 
else
    dt1 = 0.001;
end

[psi_0,mu] = get_ground_state(psi_0,dt1,g,x,k,1,V_i,N); % actually wave-function
[E_k0,E_int0,E_V0,A,B,mu0,E_tot0] = get_energy(x,psi_0,g,V_i);

%% g = 0
[T,y] = ode45(@(t,y)Ermakov1(y,1,0,wf),[0:dt:tf],[1;0]);
aa= y(:,1);
P = pi/wf;amp = 1*(wi^2-wf^2)/wf^2;
aa_ = sqrt(1+(amp).*sin(pi.*T./P).^2);
%% g = inf
[T1,y1] = ode45(@(t,y)Ermakov1(y,0,1,wf),[0:dt:tf],[1;0]);
aa_1= y1(:,1);

%P1 = exp((0.029*wf^3-0.07*wf^2+0.17*wf+1)*log(P));
P1 = exp((- 0.031*wf^4+0.098*wf^3 - 0.12*wf^2 +0.18*wf+1)*log(P));
 
amp1 =1*2*(wi^2-wf^2)/(wf^2); %pridict
aa1 = sqrt(2+amp1.*sin(pi.*T1./P1).^2)-(sqrt(2)-1);

%%
[T2,y2] = ode45(@(t,y)Ermakov1(y,A,B,wf),[0:dt:tf],[1;0]);
aa_2= y2(:,1);
aa2 = (A*aa+B*aa1);
%%
figure(1)
subplot(3,1,2)
plot(T2,aa_2,'k-','Linewidth',1.0)
hold on
plot(T,aa2,'k:','Linewidth',1.5)
hold on
set(gca,'LineWidth',1.1,'FontSize',30,'Fontname','Times New Roman');
xlabel('$t$','interpret','latex')
ylabel('$a(t)$','interpret','latex')


end

for i = 1:1:5
i
wf = 0.5;
if wf< 0.06
    tf = 200;
else
    tf = 50;
end

%% g= 1

G = [0.01,0.5,1,2,100];
g = G(i)
w_i = 1;V_i = 0.5*w_i^2.*x.^2; N=1;
a_i = sqrt(1/w_i);% harmonic length
psi_0=sqrt(N/a_i)*(1/pi)^(1/4)*exp(-x.^2./(2*a_i^2)); % Initial ansatz wavefunction

if g > 90
    dt1 = 0.0001; 
else
    dt1 = 0.001;
end

[psi_0,mu] = get_ground_state(psi_0,dt1,g,x,k,1,V_i,N); % actually wave-function
[E_k0,E_int0,E_V0,A,B,mu0,E_tot0] = get_energy(x,psi_0,g,V_i);

%% g = 0
[T,y] = ode45(@(t,y)Ermakov1(y,1,0,wf),[0:dt:tf],[1;0]);
aa= y(:,1);
P = pi/wf;amp = 1*(wi^2-wf^2)/wf^2;
aa_ = sqrt(1+(amp).*sin(pi.*T./P).^2);

%% g = inf
[T1,y1] = ode45(@(t,y)Ermakov1(y,0,1,wf),[0:dt:tf],[1;0]);
aa_1= y1(:,1);

%P1 = exp((0.029*wf^3-0.07*wf^2+0.17*wf+1)*log(P));
P1 = exp((- 0.031*wf^4+0.098*wf^3 - 0.12*wf^2 +0.18*wf+1)*log(P));
 
amp1 =1*2*(wi^2-wf^2)/(wf^2); %pridict
aa1 = sqrt(2+amp1.*sin(pi.*T1./P1).^2)-(sqrt(2)-1);

%% g= INT
[T2,y2] = ode45(@(t,y)Ermakov1(y,A,B,wf),[0:dt:tf],[1;0]);
aa_2= y2(:,1);

%P2 = (A*P+B*P1)/2;
%amp2 = (A*amp+B*amp1)/2; %pridict
%aa2 = sqrt(2+amp2.*sin(pi.*T2./P2).^2)-(sqrt(2)-1);
aa2 = (A*aa+B*aa1);
%%
figure(1)
subplot(3,1,3)
plot(T2,aa_2,'k-','Linewidth',1.0)
hold on
plot(T,aa2,'k:','Linewidth',1.5)
hold on
set(gca,'LineWidth',1.1,'FontSize',30,'Fontname','Times New Roman');
xlabel('$t$','interpret','latex')
ylabel('$a(t)$','interpret','latex')

end
%}