%++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% Submitted to :                                            +
% Prof. Dr.K. Srinivasan                                    +
%                                                           +
% Date - 23-Nov-2020                                        +
%                                                           +
% Submitted by -                                            +
$ Aditya Rana      (110118002)                              +
% Subhankar Biswas (110118084)                              +
% Vishal Mandarai  (110118100)                              +
%                                                           +
% In partial full-fillment of                               +
% Control-System (ICPC-21) final project                    +
%                                                           + 
% Deptartment of Instrumentation and Control Engineering    +
%                                                           +
% National Institute of Technology                          +
% Tiruchirappalli, India                                    +
%++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

%==========================================================================
% DC motor Control using lag compensation and PID
% To customize this code you need to:
% 1- Change the values of DC motor constants 
% 2- Change the zeros of lag compensation
% 3- Change the gains of PID controller
%==========================================================================

clc;
clear;
close all

% DC motor constants
J=0.02;
b=0.2;
kt=0.02;
ke=0.02;
R=2;
L=0.4;

% Transfer Function
num=[kt]; 
den=[J*L J*R+b*L b*R+ke*kt ];
disp('Open loop Transfer function without controller')
TF_DC=tf(num,den)
[numclp,denclp]=cloop(num,den,-1);
disp('Closed loop Transfer Function without controller')
tf(numclp,denclp)


% Open loop response without controller
step(num,den,0:0.1:5),grid on
title('Open loop response without controller')


% Closed loop response without controller
figure
step(numclp,denclp,0:0.1:5),grid on
title('Closed loop response without controller')


% State Spapce representation
A=[-(R/L) -(ke/L);(kt/J) -(b/J)];
B=[1/L; 0];
C=[0 1];
D=0;
disp('State Space representation:')
SYS_DC=ss(A,B,C,D)


% Checking Controlability and Observability
if det(ctrb(A,B))==0
    disp('          ----------> System is NOT Controllable <----------')
else
    disp('          ----------> System is Controllable <----------')
end
    
if det(obsv(A,C))==0
    disp('          ----------> System is NOT Observable   <----------')
else
    disp('          ----------> System is Observable   <----------')
end


% Drawing Root locus
figure
rlocus(num,den),grid on
title('Root Locus without controller')


% Design Criteria
Ts=1;       % Settling time<1 second

PO=0.05;    % Overshoot<5%

SSE=0.4;    % Steady state error<0.4%

abs(roots([1+(((-log(PO))/pi)^2) 0 -(((-log(PO))/pi)^2)])); % Damping ratio

Damp=ans(1); 

Wn=4/(Ts*Damp);     % Natural frequency

disp('Desired Damping ratio is:'),Damp

disp('Desired Natural Frequency is:'),Wn


% Desired Characteristic Equation:
dend=[1 2*Wn*Damp Wn^2];
disp('Desired Characteristic Equation is:'),dend


% Desired Poles location
Dp=roots(dend);
disp('Desired Pole locations:'),Dp


% From root locus and the location of desired closed loop pole, it can be
% found that a lag compensator is needed to shift the current root locus to right.
% Designing Lag compensator to meet the desired Settling time and Overshoot
% ------------------------------------------------------------------------%
z1=14;          % Assuming zero of the first lag compensator


% Finding pole of the first lag compensator
num=num/den(1);  
den=den/den(1);

ANS=inv([den(1) -dend(1) 0;den(2) -dend(2) num(1);den(3) -dend(3) num(1)*z1])*[dend(2)-den(2);dend(3)-den(3);0];

disp('Pole of the first lag compensator is:')
p1=ANS(1)
c=ANS(2);
disp('Gain of the first lag compensator is:')
K=ANS(3)


% TF of the first lag compensator G1(s)=K(s+z1)/(s+p1)
numlag1=K*[1 z1];
denlag1=[1 p1];
disp('Transfer function of the first Lag compensator to improve Ts and PO%:')
tf(numlag1,denlag1)


% DC motor Transfer function with Lag compensator 1
disp('DC motor Transfer function with Lag compensator')
NUM=conv(numlag1,num);
DEN=conv(denlag1,den);
TF=tf(NUM,DEN)


% Root locus with Lag compensator 1
figure
rlocus(TF),grid on
title('Root locus with Lag compensator 1')


% Open loop responce of the system with Lag compensator 1
figure
step(TF,0:0.1:5),grid on
title('Open loop response with lag compensator 1')


% Closed loop responce of the system with Lag compensator 1
[numc,denc]=cloop(NUM,DEN);
figure
step(numc,denc,0:0.1:5),grid on
title('Closed loop response with Lag compensator 1 that improves Ts & PO%')


% Improving SSE by adding a second lag compensator
z2=2.9;     % Assuming zero of the 2nd lag compensator
SSE=0.004;  % Steady State Error design criteria


% Solving for pole of the 2nd lag compensator
disp('pole of the 2nd lag compensator')
p2=(1+((K*z1*num(1)/denlag1(2))/den(3)))*z2*SSE
numlag2=[1 z2];
denlag2=[1 p2];
NumLag=conv(numlag1,numlag2);
DenLag=conv(denlag1,denlag2);
disp('The 2nd Lag compensator Transfer function to improve SSE:')
tf(numlag2,denlag2)
disp('The overal Lag compensator transfer function (lag1*lag2):')
tf(NumLag,DenLag)


% DC motor transfer function with Lag compensator that improves Ts, PO% & SSE
NumDC=conv(NumLag,num);
DenDC=conv(DenLag,den);
disp('Open loop TF of the DC motor with final Lag compensator (improved Ts, PO% & SSE) ')
tf(NumDC,DenDC)


% Root locus with final lag compensator
figure
rlocus(NumDC,DenDC), grid on
title('Root locus with final lag compensator')


% Closed loop TF of the DC motor with Lag compensator
[NumCLP,DenCLP]=cloop(NumDC,DenDC);
disp('closed loop TF of the DC motor with final Lag compensator (improved Ts, PO% & SSE) ')
tf(NumCLP,DenCLP)
figure
step(NumCLP,DenCLP,0:0.1:5), grid on
title('Closed loop response with final Lag compensator')


%--------------------End of Lag compensator Design------------------------%
%--------------------------------PID control------------------------------%
% PID control gain, using trial and error
kp=70;
ki=170;
kd=5;


% PID transfer function
numPID=[kd kp ki];
denPID=[1 0];


% Open loop TF of DC motor with PID controller
num_DC_PID=conv(num,numPID);
den_DC_PID=conv(den,denPID);
disp('Open loop TF of DC motor with PID controller')
tf(num_DC_PID,den_DC_PID)


% Closed loop TF of DC motor with PID controller
[NumPID_CLP,DenPID_CLP]=cloop(num_DC_PID,den_DC_PID);
disp('Closed loop TF of DC motor with PID controller')
tf(NumPID_CLP,DenPID_CLP)
figure
step(NumPID_CLP,DenPID_CLP), grid on
title('Closed loop response of DC with PID Control')


%-------------------------End of PID control------------------------------%
% Bode plot, Determining gain and phase margin
figure
margin(numclp,denclp), grid on

figure
margin(numc,denc), grid on %Bode plot of closed loop TF with lag compensator

figure
margin(NumPID_CLP,DenPID_CLP), grid on % Bode plot of closed loop TF with PID controller