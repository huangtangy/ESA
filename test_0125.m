clc;close 
clf;clear
syms t a
hbar=1;m=1;%Dimensionless and defined the coordinate and kintic space
n=2^15;dx=0.01;x=(-n/2:n/2-1)*dx;
dk=2*pi/(n*dx); M =n/2;Nx = 2*M;
k =(-n/2:n/2-1)*dk;%Define k-space grid
N = 1;
tic

%%
shell_A1 = zeros();shell_B1 = zeros();shell_g1 = zeros();
shell_tf = zeros();shell_a1 = zeros();shell_a2 = zeros();
shell_psi0 = zeros();shell_psi1=zeros();shell_psi2=zeros();
EE = zeros();shell_Bures = zeros();shell_QSL = zeros();
shell_control = zeros();
shell_fidelity=zeros();
dt = 0.001;%time grid
Ni = 50;


AA=[0.999341056441508,0.999243267952929,0.999131006214376,0.999002127273517,0.998854179467035,0.998684339927594,0.998489379164041,0.998265585811989,0.998008703459152,0.997713854930085,0.997375439799454,0.996987044366654,0.996541310871890,0.996029810237536,0.995442887387883,0.994769478809948,0.993996918921688,0.993110717169562,0.992094293584866,0.990928693158820,0.989592264245005,0.988060283638185,0.986304552381463,0.984292936030572,0.981988856543039,0.979350736602220,0.976331381589263,0.972877320328731,0.968928097618334,0.964415515598615,0.959262867767580,0.953384141053749,0.946683273451462,0.939053438282240,0.930376476240541,0.920522509930813,0.909349837710588,0.896705259144915,0.882424942755549,0.866336052011023,0.848259341578257,0.828012936480639,0.805417620455340,0.780303826226824,0.752520622256885,0.721946813751401,0.688504179225486,0.652172651268161,0.613006864732884,0.571153148014631,0.526655402956935,0.480296905390318,0.432383067077457,0.383543328095823,0.334521595993201,0.286151021433003,0.239315091319584,0.194896639785719,0.153718536712543,0.116482643496494,0.0837150022163788,0.0557256179941182,0.0325896481452531,0.0141535053583162,6.51126601086379e-05,-0.0101768507335671,-0.0171631211532643,-0.0215131575601821,-0.0238223263376342,-0.0246232976864755,-0.0243637690875665,-0.0233998528952066,-0.0220012260216013,-0.0203634654541645,-0.0186231189860438,-0.0168722988353843,-0.0151710044167637,-0.0135566265474863,-0.0120508705987214,-0.0106646526203871,-0.00940151970118852,-0.00826005617164380,-0.00723558663744113,-0.00632140277187818,-0.00550963193943121,-0.00479185094464997,-0.00415951761424052,-0.00360427012302930,-0.00311811385268551,-0.00269354625852556,-0.00232362091850126,-0.00200197904310632,-0.00172284771885687,-0.00148102759544539,-0.00127186352924769,-0.00109120887647393,-0.000935388040332178,-0.000801154548764178,-0.000685651200811417,-0.000586371992398739,-0.000586371992398739];
BB=[0.000660746880298286,0.000758598077460778,0.000870933288546918,0.000999894239643011,0.00114793864693848,0.00131788661352555,0.00151297376957191,0.00173691210883836,0.00199395960622028,0.00228899985465084,0.00262763309102342,0.00301628019193361,0.00346230138252314,0.00397413162185854,0.00456143486663421,0.00523527959323248,0.00600833828898686,0.00689511378013822,0.00791219552731558,0.00907854923204479,0.0104158432679299,0.0119488154089257,0.0137056834766880,0.0157186031025017,0.0180241754401354,0.0206640068849273,0.0236853213499101,0.0271416239998528,0.0310934123770837,0.0356089272228399,0.0407649297679595,0.0466474857100203,0.0533527267550219,0.0609875488308899,0.0696701914600452,0.0795306245051101,0.0907106457597953,0.103363567449743,0.117653339652288,0.133752926139075,0.151841716416814,0.172101728535458,0.194712339211177,0.219843275390780,0.247645635921147,0.278240786114579,0.311707111568008,0.348064846083020,0.387259508865167,0.429144919009665,0.473377412072219,0.519739572742589,0.567657393820690,0.616501463322473,0.665527889177756,0.713903564422674,0.760744967841948,0.805169349674488,0.846353872134980,0.883596718931271,0.916371852783235,0.944369401460319,0.967514233660793,0.985960023504268,1.00005895656930,1.01031244454427,1.01731131119598,1.02167515754156,1.02399940640938,1.02481693494824,1.02457554081043,1.02363149619478,1.02225465474690,1.02064078849930,1.01892664482494,1.01720455529662,1.01553476834713,1.01395494669308,1.01248708661567,1.01114243082337,1.00992488355180,1.00883342234787,1.00786379319991,1.00700975780879,1.00626396126229,1.00561854443284,1.00506558077960,1.00459739093425,1.00420672283296,1.00388689131983,1.00363184541738,1.00343621136250,1.00329529301887,1.00320507614406,1.00316220523157,1.00316395960834,1.00320823101153,1.00329349170313,1.00341877217451,1.00358364139575,1.00358364139575];
% the actual width
shell_a11 = [1.00009961168818,1.00011437795046,1.00013133171810,1.00015079674086,1.00017314528591,1.00019880406721,1.00022826354999,1.00026208644657,1.00030091900592,1.00034550281461,1.00039668935991,1.00045545625734,1.00052292541648,1.00060038463923,1.00068931231136,1.00079140522830,1.00090861062963,1.00104316333995,1.00119762829030,1.00137494812261,1.00157849947737,1.00181215641085,1.00208036257579,1.00238821511481,1.00274155947984,1.00314709778282,1.00361251322580,1.00414661160320,1.00475948271775,1.00546268339900,1.00626944712571,1.00719491973225,1.00825642823392,1.00947378402018,1.01086962416072,1.01246979540353,1.01430378362289,1.01640519292141,1.01881227576092,1.02156851663334,1.02472327020156,1.02833245041295,1.03245926473398,1.03717498635152,1.04255974681923,1.04870332450331,1.05570590140397,1.06367873995731,1.07274472551167,1.08303870759877,1.09476203772598,1.10797294588062,1.12288821168833,1.13968718477871,1.15855696164441,1.17968991064040,1.20328068146568,1.22952294851437,1.25860579022556,1.29071068024248,1.32600877604066,1.36465914636339,1.40680814901226,1.45258989817576,1.50212778909822,1.55553681114243,1.61292646052461,1.67440356593066,1.74007524388204,1.81005144296966,1.88444725206540,1.96338464101821,2.04699403136458,2.13541562805800,2.22880039239829,2.32731089552695,2.43112196383986,2.54042137846433,2.65541039004487,2.77630418560068,2.90333234592513,3.03673929119815,3.17678470765424,3.32374401486077,3.47790886988033,3.63958751828110,3.80910541910838,3.98680578539329,4.17305007634300,4.36821866764795,4.57271148519115,4.78694868611359,5.01137139673113,5.24644239799668,5.49264697890545,5.75049375642129,6.02051564669238,6.30327055006358,6.59934243041495,6.90934225592258];
shell_a22 = [3.16327465825084,3.16342237402724,3.16359195796511,3.16378664297074,3.16401017391391,3.16426678135913,3.16456139650184,3.16489960920823,3.16528788857068,3.16573364718475,3.16624535167020,3.16683277594158,3.16750708503490,3.16828113194951,3.16916962444743,3.17018943735492,3.17135991530457,3.17270325883052,3.17424493628819,3.17601411733342,3.17804412086202,3.18037325599135,3.18304531666921,3.18611041585117,3.18962596949364,3.19365736609300,3.19827958116302,3.20357802857262,3.20965020036929,3.21660714070077,3.22457512793034,3.23369773959763,3.24413790181208,3.25607987914991,3.26973195991471,3.28532871433158,3.30313366186581,3.32344161262441,3.34658109198739,3.37291636280837,3.40284895033487,3.43681835315777,3.47530203102632,3.51881386466314,3.56790132517726,3.62314040886256,3.68512931567285,3.75447993146240,3.83180786822895,3.91772194491634,4.01285103622956,4.11769147840231,4.23280873093086,4.35868273479850,4.49576412907344,4.64447552234143,4.80518604362896,4.97826199611995,5.16404163132313,5.36285944547264,5.57506517260771,5.80100169156697,6.04105638010954,6.29560849332006,6.56508389739221,6.84995126325218,7.15068543333744,7.46782201560809,7.80193460474431,8.15362143788827,8.52353708321080,8.91238929317673,9.32091068728153,9.74990902670626,10.2002065051553,10.6727143122065,11.1683683938989,11.6881702728533,12.2331750184265,12.8044922486289,13.4032985429232,14.0308182684336,14.6883726903604,15.3773052422604,16.0990612515563,16.8551332173149,17.6471202665163,18.4768382824960,19.3502090420633,20.2878552751483,21.1983592575423,21.9857783350397,22.6543195257572,23.2224770058279,23.1201316453345,24.1211408114988,24.4764128774473,21.8634753312044,21.4370027157661,21.2690780495265];
shell_A11 = [0.999341055835202,0.999243268798177,0.999131004895832,0.999002128461845,0.998854178940196,0.998684340213325,0.998489379827878,0.998265585642154,0.998008703440546,0.997713853735199,0.997375440033039,0.996987044527119,0.996541310909037,0.996029811531701,0.995442887784119,0.994769478493176,0.993996919605751,0.993110716749321,0.992094292406391,0.990928692572580,0.989592264400601,0.988060283430993,0.986304552674759,0.984292936431335,0.981988856525065,0.979350735690765,0.976331380548717,0.972877320096235,0.968928097610494,0.964415517420200,0.959262868136578,0.953384141449503,0.946683273667649,0.939053438313644,0.930376477551712,0.920522510508018,0.909349837824101,0.896705259486737,0.882424942314708,0.866336052654618,0.848259341605992,0.828012936599981,0.805417619922556,0.780303826992928,0.752520622273118,0.721946812975100,0.688504180121192,0.652172651207296,0.613006864177981,0.571153148236772,0.526655421168361,0.480296919327768,0.432383077312730,0.383543332508030,0.334521597012890,0.286151019904735,0.239315105644054,0.194896641347100,0.153718536852615,0.116482642935415,0.0837150022233910,0.0557256182437846,0.0325896457105128,0.0141535025618622,6.51123703196992e-05,-0.0101768533849073,-0.0171631213447073,-0.0215131489816272,-0.0238223285625516,-0.0246232990687924,-0.0243637678407575,-0.0233998492743510,-0.0220012261700779,-0.0203634676130350,-0.0186231211563804,-0.0168722987364146,-0.0151710037965553,-0.0135566258579951,-0.0120508704305019,-0.0106646533620820,-0.00940152059097482,-0.00826005514951420,-0.00723558565812902,-0.00632140350157152,-0.00550963220397103,-0.00479185065862477,-0.00415951821788111,-0.00360427020880619,-0.00311811387894129,-0.00269354621974574,-0.00232362131839087,-0.00200197911582120,-0.00172284757255670,-0.00148102744460268,-0.00127186342366944,-0.00109120933617632,-0.000935388340582931,-0.000801154535912832,-0.000685651246362455,-0.000586372215817955];
shell_B11 = [0.000660746880314149,0.000758598077394548,0.000870933288470076,0.000999894239926589,0.00114793864682920,0.00131788661373872,0.00151297376985158,0.00173691210876820,0.00199395960627597,0.00228899985470065,0.00262763309165419,0.00301628019209783,0.00346230138130691,0.00397413162235235,0.00456143486661454,0.00523527959262994,0.00600833828991327,0.00689511378245746,0.00791219552496052,0.00907854923199679,0.0104158432685725,0.0119488154091179,0.0137056834787486,0.0157186031023145,0.0180241754406249,0.0206640068827330,0.0236853213449940,0.0271416239932707,0.0310934123774540,0.0356089272409659,0.0407649297815570,0.0466474857236818,0.0533527267693424,0.0609875488279431,0.0696701914660107,0.0795306245060732,0.0907106457712799,0.103363567488871,0.117653339662050,0.133752926125817,0.151841716410067,0.172101728555902,0.194712339197191,0.219843275460582,0.247645635914405,0.278240786032078,0.311707111595074,0.348064846004714,0.387259508860945,0.429144919128267,0.473377424215703,0.519739584507429,0.567657404311856,0.616501467340741,0.665527891034724,0.713903560250920,0.760745014255573,0.805169359377148,0.846353869283437,0.883596704328831,0.916371857698947,0.944369404215131,0.967514246935706,0.985960042168264,1.00005895816108,1.01031245394650,1.01731131315493,1.02167512913307,1.02399941239450,1.02481694103436,1.02457553719421,1.02363148504582,1.02225465503928,1.02064079483836,1.01892665109951,1.01720455529662,1.01553476664989,1.01395494400258,1.01248708650045,1.01114243304760,1.00992488599068,1.00883341906933,1.00786379011087,1.00700975985470,1.00626396213645,1.00561854359208,1.00506558263932,1.00459739136272,1.00420672289331,1.00388689121902,1.00363184647503,1.00343621170742,1.00329529277757,1.00320507552487,1.00316220469346,1.00316396103920,1.00320823172153,1.00329349157301,1.00341877243421,1.00358364208671];




for i = 1:1:Ni
   
i
if i>Ni/2
    dt = 0.001;
end

j = -3+(i)*(6/Ni);g =10^(j)*0.5;
w_i = 1;a_i = (1/w_i)^(1/2);
w_f = 0.1;a_f = (1/w_f)^(1/2);

%% producing the ground states
A = AA(i);B=BB(i);
a_i=max(double(vpasolve(a*w_i^2-A/(a^3)-B/(a^2)==0,a,[0,inf])));%exact initial width
a_f=max(double(vpasolve(a*w_f^2-A/(a^3)-B/(a^2)==0,a,[0,inf])));%exact final width

g_i = g
V_i = 0.5*m*w_i^2.*x.^2; %initial parameters frequency and interaction
psi_i=sqrt(N/a_i)*(1/pi)^(1/4)*exp(-x.^2/(2*a_i^2)); %Initial ansatz wavefunction

[psi_0,mu] = get_ground_state(psi_i,dt,g_i,x,m,V_i,N); % actually wave-function
shell_psi0(i,1:length(psi_0)) = psi_0;

psi_f=sqrt(N/a_f)*(1/pi)^(1/4)*exp(-x.^2/(2*a_f^2)); %final ansatz wavefunction
V_f = 0.5*w_f^2.*x.^2; 
[psi_1,muf] = get_ground_state(psi_f,dt,g,x,1,V_f,N); % actually wave-function
shell_psi1(i,1:length(psi_1)) = psi_1;

%% calculating A and B
%{
n0 = (abs(psi_0).^2); % density
sig_s = trapz(x,x.^2.*n0); % sigma squre
E_int = 0.5*g_i*trapz(x,n0.^2); 
x1=((-M+1):1:M-1).*dx; 
E_k = 0.5*trapz(x1,(diff(sqrt(n0))./dx).^2);
E_V =trapz(x,V_i.*n0);
 
dd = diff(n0,2)./dx^2;
D_0 = -(0.25./(a_i^2.*n0(M))).*dd(M);

A = ((2/(sig_s*w_i^2))*(D_0-E_k));
B = ((2/(sig_s*w_i^2))*(g_i*n0(M)-2*E_int));

%g_i*n0(M)+D_0;%mu
%AB = ((2/(sig_s*w0^2))*(muu-E_k-2*E_int));
%delta = abs(muu-E);

shell_A1(1,i) = A; 
shell_B1(1,i)= B;
shell_g1(1,i) = g_i;
A+B
%}
%a_i=width(psi_0,x,dx);%exact initial width
%a_f=width(psi_1,x,dx);%exact final width
%}

%% computing the time for bang-bang

shell_a1(1,i) = a_i;
shell_a2(1,i) = a_f;

beta = a_i;gamma =a_f;
   
d = w_i^2;% delta
C1 = -d*beta^2 + A/(beta^2)+2*B/(beta);
C2 = d*gamma^2 + A/(gamma^2)+2*B/(gamma);%integrated cofficienr
x1B=sqrt((C2-C1)/(2*d)); %xB B

fun1= @(s)sqrt(1./(C1+d.*s.^2-A./(s.^2)-2*B./(s)));
t1 = integral(fun1,beta,x1B); % time for first segment

fun2= @(s)sqrt(1./(C2-d.*s.^2-A./(s.^2)-2*B./(s)));
t2 = integral(fun2,x1B,gamma); % time for second segment

tf = t1+t2; %total time
shell_tf(1,i) = tf;


%% time evolution for IE
%{
dt=10^(-2); 
tf3= 1;
Nt = tf3/dt;
T = 0:tf3/Nt:tf3;% total time can be arbitary 
a1(t)=beta-6*(beta-gamma).*(t/tf3).^5+15*(beta-gamma).*(t/tf3).^4-10.*(beta-gamma).*(t/tf3).^3;%polynomial 
a2 = diff(a1,t,2); % ddot(a)
w2 = A./a1(T).^4+B./(a1(T).^3)-a2(T)./a1(T);%omega^2
psi = psi_0;
for itime=0:Nt-1 %Time-stepping with split-step Fourier method 
    
     u =double(w2(itime+1));
     V=0.5*m*u*x.^2/hbar; %Define potential
     psi = FFT( psi,V,g,dt,k,x);
     %plot(x,abs(psi).^2,x,abs(psi_0).^2,x,abs(psi_T).^2)
     %drawnow
end

%}
%% time evolution for bang-bang control

dt = 0.01;
psi = psi_0;Nt = round(tf/dt);
d = w_i^2;
for itime=0:Nt %Time-stepping with split-step Fourier method 
    t = itime*dt;
    
    if t ==0
        u = w_i^2;
    elseif (t>0)&&(t<t1)
        u = -d;
    elseif (t>=t1)&&(t<tf)
        u = d;
    elseif t ==tf 
        u = w_f^2;
    end
    shell_control(1,itime+1) = u;
    V=0.5*m*u*x.^2/hbar; %Define potential
    psi = FFT( psi,V,g,dt,k,x); 
    
end
%}
%%
shell_psi2(i,1:length(psi)) = psi;
fidelity = abs((sum((conj(psi).*psi_1)).*dx)).^2/N
shell_fidelity(1,i) = fidelity;
Bures = acos(abs((sum(conj(psi).*psi_0).*dx)));
shell_Bures(1,i) = Bures;

%% trajectories and time-averagy energy, QSL
%first seagma
delta = -1;
[T1,y1] = ode45(@(t,y1)Ermakov1(y1,A,B,delta),[0:dt:t1],[a_i;0]);
a1= y1(:,1);
da1 = y1(:,2);
E1 = 0.5.*da1+A./(2.*a1.^2)+B./a1+0.5.*(delta).*a1.^2;

% second seagma
delta = 1;
[T2,y2] = ode45(@(t,y2)Ermakov1(y2,A,B,delta),[0:dt:t2],[a1(end);da1(end)]);
a2= y2(:,1);
da2 = y2(:,2);
E2 = 0.5.*da2+A./(2.*a2.^2)+B./a2+0.5.*(delta).*a2.^2;

E_ave = (trapz(T1,E1)+trapz(T2,E2))./tf;
qsl = (sin(Bures))^2/(2*E_ave);
shell_QSL(1,i) = qsl;

end

plot(x,abs(psi).^2,x,abs(psi_1).^2)%,x,abs(psi_0).^2)
%plot([0:dt:tf],shell_control)
%plot(T1,a1,T2,a2);

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

G = linspace(0.001,1000,length(AA)-1);
%{
subplot(2,1,1)
plot(G,shell_tf,'Linewidth',1.0)
xlabel('g')
ylabel('$t_f$','interpret','latex')
hold on
set(gca,'LineWidth',1.1,'FontSize',20,'Fontname','Times New Roman');

subplot(2,1,2)
plot(G,shell_a1,G,shell_a2,'Linewidth',1.0)
%plot(G,shell_a11,'--',G,shell_a22,'--','Linewidth',1.0)
xlabel('g')
ylabel('$a_1/a_2$','interpret','latex')
%plot(shell_g1,shell_A1+shell_B1,'-.k','Linewidth',1.0)
set(gca,'LineWidth',1.1,'FontSize',20,'Fontname','Times New Roman');
%}
toc