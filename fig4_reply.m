clc;close 
clf;clear
syms a t
hbar=1;m=1;%Dimensionless
dx=double(1e-2);M=5000;x=(-M:1:M).*dx;Nx=2*M+1;%step split
dk=pi/(M*dx); k=(-M:1:M)*dk; N =1;%Define k-space grid
w_0 = 250*2*pi;w_i = 250*2*pi/w_0;w_f = 1/10^2;%intial and final frequency 
F1=zeros();TTime1=zeros();G1=zeros(); Z1 = zeros();
F1_1 =zeros();TTime1_1=zeros();G1_1=zeros();
F2=zeros();TTime2=zeros();G2=zeros(); Z2 = zeros();
F2_1=zeros();TTime2_1=zeros();G2_1=zeros(); Z2_1 = zeros();
F3=zeros();TTime3=zeros();G3=zeros(); Z3 = zeros();
F3_1=zeros();TTime3_1=zeros();G3_3=zeros();
tic
%%  %-------------------------------bang-medthod------------------------------%
w_i = 1;w_f = 0.2;
%jj = 1;
%{
for jj =1:50
    gN = 0+10/50*(jj-1)
    jj
    %a_i = 1;a_f =10;
    a_i=double(vpasolve(a*w_i^2-1/(a^3)==0,a,[0,inf]));%exact initial width #-gN/(a^2*sqrt(pi*2))
    a_f=double(vpasolve(a*w_f^2-1/(a^3)==0,a,[0,inf]));%exact final width
    beta = a_i;gamma =a_f;
    psi_0=sqrt(N/a_i)*(1/pi)^(1/4)*exp(-x.^2/(2*a_i^2)); %Initial ansatz wavefunction
    psi_T=sqrt(N/a_f)*(1/pi)^(1/4)*exp(-x.^2/(2*a_f^2)); %final ansatz wavefunction 

    w1 = sqrt(1/(beta^2*gamma^2)+2*gN/(sqrt(2*pi)*beta*gamma*(beta+gamma)));%w_c inter omega
    U_i = w1^2*beta^2/2+1/(2*beta^2)+gN/(sqrt(2*pi)*beta);% initial effective potential
    fun = @(a)1./sqrt(2*(U_i-(w1^2*a.^2./2+1./(2*a.^2)+gN./(sqrt(2*pi).*a))));
    tf1 = integral(fun,beta,gamma);  % time comsuer
    dt = 10^(-2);Nt = tf1/dt;
    %Nt=100;dt=tf1/Nt;  %psi= psi_0;%time split 
    %[psi_0,u1] = get_ground_state(psi_0,dt,gN,x,k,m,0.5*w_i^2*x.^2,N);%exact initial wavefunction
    %V=0.5*m*w_i^2*x.^2/hbar;N0=1;
    %[psi_T,u2] = get_ground_state(psi_T,dt,gN,x,k,m,0.5*m*w_f^2*x.^2/hbar,N);%exact final wavefunction
    psi = psi_0;
for itime=0:Nt %Time-stepping with split-step Fourier method 
    t = itime*dt;
    if t ==0
        u=w_i^2;
    elseif (t>0)&&(t<tf1)
        u=w1^2;
    elseif t ==tf1
        u=w_f^2;
    end
     V=0.5*m*u*x.^2/hbar; %Define potential
     psi = FFT( psi,V,gN,dt,k,x); 
     %plot(x,abs(psi).^2,x,abs(psi_0).^2,x,abs(psi_T).^2)
     %drawnow
end
TTime1(1,jj)=tf1; 
G1(1,jj) = gN;
F1(1,jj) = abs((sum(conj(psi).*psi_T).*dx)).^2/N

end
%}
%axis([0 1.01 0.9 10.3])
%text(0.95,9,{'(b)'})

%%  %-------------------------------bangbang-medthod---------------------------%
%{
PP = zeros();Time = zeros();
jj = 1;
for jj = 1:50 % gN=[-1,1]
    gN = -1+2/50*jj; %4/100*jj-2
    %a_i = 1;a_f =10;
    a_i=double(vpasolve(a*w_i^2-1/(a^3)+gN/(a^2*sqrt(pi*2))==0,a,[0,inf]));%exact initial width
    a_f=double(vpasolve(a*w_f^2-1/(a^3)+gN/(a^2*sqrt(pi*2))==0,a,[0,inf]));%exact final width
    beta = a_i;gamma =a_f;
    psi_0=sqrt(N/a_i)*(1/pi)^(1/4)*exp(-x.^2/(2*a_i^2)); %Initial ansatz wavefunction
    psi_T=sqrt(N/a_f)*(1/pi)^(1/4)*exp(-x.^2/(2*a_f^2)); %final ansatz wavefunction 
    
    d = w_i^2;% delta
    C1 = -d*beta^2 + 1/(beta^2)+2*gN/(sqrt(2*pi)*beta);
    C2 = d*gamma^2 + 1/(gamma^2)+2*gN/(sqrt(2*pi)*gamma);%integrated cofficienr
    x1B=sqrt((C2-C1)/(2*d)); %xB B
    fun1= @(s)sqrt(1./(C1+d.*s.^2-1./(s.^2)-2*gN./(sqrt(2*pi).*s)));
    t1 = integral(fun1,beta,x1B); % time for first segment
    fun2= @(s)sqrt(1./(C2-d.*s.^2-1./(s.^2)-2*gN./(sqrt(2*pi).*s)));
    t2 = integral(fun2,x1B,gamma); % time for second segment
    tf2 = t1+t2; %total time
    Nt = 1000;dt = tf2/Nt;
    %Nt=1000;dt=tf2/Nt;  
    %[psi_0,u1] = get_ground_state(psi_0,dt,gN,x,k,m,0.5*w_i^2*x.^2,N);%exact initial wavefunction
    %[psi_T,u2] = get_ground_state(psi_T,dt,gN,x,k,m,0.5*w_f^2*x.^2,N);%exact final wavefunction
    psii = psi_0;
    ii = 1;
for itime=0:Nt %Time-stepping with split-step Fourier method 
    t = itime*dt;
    if t ==0
        u=w_i^2;
    elseif (t>0)&&(t<=t1)
        u=-d;
    elseif (t>t1)&&(t<tf2)
        u=d;
    elseif t == tf2
        u=w_f^2;%w_f^2;
    end
     V=0.5*m*u*x.^2/hbar; %Define potential
     psi = FFT(psii,V,gN,dt,dx,M);
     ss = abs(psi).^2;
     
     %for i= 1:length(psi)
         %PP(ii,i) = ss(i);
     %end
     
     
     psii = psi;
     %plot(x,abs(psi).^2,x,abs(psi_T).^2)%,x,abs(psi_0).^2
     %drawnow
     Time(ii,1) = t;
     ii = ii+1;
     
end
TTime2(1,jj)=tf2; 
G2(1,jj) = gN;
F2(1,jj) = abs((sum(conj(psi).*psi_T).*dx)).^2/N
%jj = jj +1;
end

figure
[X,Y] = meshgrid(x,Time./tf2);
mesh(X,Y,PP)
set(gca,'LineWidth',1.1,'FontSize',22,'Fontname','Times New Roman');
ylabel('t/t_f')
xlabel('x')
zlabel('|\psi|^2')
%}
%%  %-----------------------------------IE-medthod-----------------------------%

%jj = 1;
%a_i = 1;a_f =10;

for jj = 1:10 % gN=[-1,1]
    gN =0+10/20*(jj-1)
    %a_i =1; a_f =10;
    a_i=double(vpasolve(a*w_i^2-1/(a^3)==0,a,[0,inf]));%exact initial width -gN/(a^2*sqrt(pi*2))
    a_f=double(vpasolve(a*w_f^2-1/(a^3)==0,a,[0,inf]));%exact final width
    beta = a_i;gamma =a_f;
    psi_0=sqrt(N/a_i)*(1/pi)^(1/4)*exp(-x.^2/(2*a_i^2)); %Initial ansatz wavefunction
    psi_T=sqrt(N/a_f)*(1/pi)^(1/4)*exp(-x.^2/(2*a_f^2)); %final ansatz wavefunction 
    
    
    dt=10^(-2); 
   %Nt=5000; %step
    tf3= 1;
    %[tf3 ww] =TimeforIE(6.5,gN,500);aa=ww;
    
    Nt = tf3/dt;
    T = 0:tf3/Nt:tf3;% total time can be arbitary 
    a1(t)=beta-6*(beta-gamma).*(t/tf3).^5+15*(beta-gamma).*(t/tf3).^4-10.*(beta-gamma).*(t/tf3).^3;%polynomial 
    a2 = diff(a1,t,2); % ddot(a)
    w2 = 1./a1(T).^4+gN./(sqrt(2*pi).*a1(T).^3)-a2(T)./a1(T);%omega^2

     %psi= psi_0;%time split
    [psi_0,u1] = get_ground_state(psi_0,dt,gN,x,k,m,0.5*w_i^2*x.^2,N);
    [psi_T,u2] = get_ground_state(psi_T,dt,gN,x,k,m,0.5*w_f^2*x.^2,N);
    psi = psi_0;
for itime=0:Nt-1 %Time-stepping with split-step Fourier method 
    
     u =double(w2(itime+1));
     V=0.5*m*u*x.^2/hbar; %Define potential
     psi = FFT( psi,V,gN,dt,k,x); 
     %plot(x,abs(psi).^2,x,abs(psi_0).^2,x,abs(psi_T).^2)
     %drawnow
end

TTime3(1,jj)=tf3; 
G3(1,jj) = gN;
F3(1,jj) = abs((sum(conj(psi).*psi_T).*dx)).^2/N
    %jj = jj +1;
end
%}

plot(G1,F1,'k-')
plot(G3,F3,'k--')
xlabel('$gN$','interpret','latex')
ylabel('$Fidelity$','interpret','latex')


%%% omega_f = 0.2
%%bangbang
%%IE





%}