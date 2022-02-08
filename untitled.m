
%% checking the anlytical expression of quench dynamic
wf = 0.1;
tf = pi/(2*wf);

dt = 0.01;
T = 0:dt:tf;
A = 1;B =0;
[T,y] = ode45(@(t,y)Ermakov1(y,1,0,wf),[0:dt:tf],[1;0]);
at = y(:,1);
At = sqrt(1+(1-wf^2).*(sin(wf.*T)).^2/(wf^2));

plot(T,At,'r-',T,at,'b--')