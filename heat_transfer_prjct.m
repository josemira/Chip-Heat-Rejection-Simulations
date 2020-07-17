%% Inputs
clear
clc
 
Tc = 85; %Celsius

Tinf = 20; %Celsius

h = 100; %W/m^2 K
k = 180; %W/m K

Lf = .012; % meter
Lb = .003; % meter

L = Lf;

W = 20;%in millimeters
N=3:1:10;
t = ((W - 2.*N + 2)./N)./1000;

%re-assigning W value
W = 0.02;%in meters

P = 2.*(W+t);
Ac = t*W;
Af = 2.*((Lf.*W)+(Lf.*t));
At = W.^2 - N.*(W.*t) + N.*Af;
%% N vs qc

figure(1)
m = sqrt((h.*P)./(k.*Ac));

nf = (tanh(m.*L))./(m.*L);
no = 1-((N.*(Af./At)).*(1-nf));

Rs = (2*10^-5);
Rc = Lb/k;
Rf = 1./(no.*h.*At);

qc = (Tc-Tinf)./(Rs+Rc+Rf);
plot(N,qc)
xlabel('number of fins')
ylabel('heat flux (W/m^2)')
title('N vs qc')
%% N,Tb
figure(2)
Tb = qc.*Rf + Tinf;
plot(N,Tb)
xlabel('number of fins')
ylabel('Base Temp "Tb" (C)')
title('N vs Tb')
%% N,qf
figure(3)

theta_b = Tb - Tinf;
M = sqrt(h.*P.*k.*Ac.*theta_b);
qf = M.*tanh(m.*L);
%qf = nf.*h.*Af.*(Tb-Tinf);
plot(N,qf)
xlabel('number of fins')
ylabel('heat flux from fins (W/m^2)')
title('N vs qf')
%% N,Ttip
figure(4)
L = Lf;
x = Lf;
Ttip = Tinf + (Tb-Tinf).*cosh(m.*(L-x))./(cosh(m.*L));
plot(N,Ttip)
xlabel('number of fins')
ylabel('Tip Temp "Ttip" (C)')
title('N vs Ttip')
%% wt,N
figure(5)
d = 2710;%kg/m^3 density of alluminum
Vf = N.*t.*(Lf*W); %m^3Volume of fins
mf = d.*Vf; %kg mass of fins
wt = mf.*9.81;%Newtons weight of fins

plot(N,wt)
xlabel('number of fins')
ylabel('total weight of fins (Newtons)')
title('N vs weight')
%% Adiabatic Tip Condition

%theta = T - T_inf;
%theta_b = Tb - T_inf;
%m = sqrt(h*P/k*Ac);
%M = sqrt(h*P*k*Ac*theta_b);

% Temperature Distribution
%theta_L = cosh(m*(L-x))/cosh(m*L)

%Fin Heat Transfer Rate 
%qf = M*tanh(m*L)
%% Plotting temps for fdm c/f and 1D

% T1 = Tc - qc(1,1)*Rs;
% Ttip = Ttip(1,1);

fdmCt = xlsread('fdmc.xlsx','B5:M5');
fdmFt = xlsread('fdmf.xlsx','B5:AW5');

figure(6)
hold on
%FDM Coarse temp gradient
xC = linspace(0,Lf,12);
plot(xC,fdmCt)

%FDM Fine temp gradient
xF = linspace(0,Lf,48);
plot(xF,fdmFt)

%1D temp gradient
m = m(1,1);
L = Lf;
x1D = linspace(0,Lf,11);
x = x1D;
Tb = Tb(1,1);
Tx = Tinf + (Tb-Tinf).*cosh(m.*(L-x))./cosh(m.*L);

plot(x,Tx)


title('Fin Temp gradient')
ylabel('Temp (C)')
xlabel('x position (meters)')

legend('FDM - Coarse','FDM - Fine','1D')

% -----------------------------------------------------
% q_conduction heat Transfer plots
figure(7)
hold on 

% 1D heat Transfer
% qc1d = ones(1,10);
% n=0;
% for i = 1:length(Tx)-1
%     n = n+1;
%     qc1d(n) = (k/.001).*(Tx(1,n) - Tx(1,n+1));
% end
% 
% qc1 = ones(1,10);
% a = 0;
% for i = 1:length(qc1d)
% a = a +1;
% qc1(a) = (0.5*.02*.001).*(qc1d(1,a) - qc1d(1,a+1));
% end
% 
onedq = xlsread('1d.xlsx','C22:L22');
xd = linspace(0,Lf,10);

plot(xd,onedq)

%FDM Coarse heat transfer
fdmCq = xlsread('fdmc.xlsx','B11:L11');
xC = linspace(0,Lf,11);
plot(xC,fdmCq)

%FDM Fine heat transfer
fdmFq = xlsread('fdmf.xlsx','B11:AV11');
xF = linspace(0,Lf,47);
plot(xF,fdmFq)  
    
title('Conduction Heat Transfer Across Fin')
ylabel('qc (W/m^2)')
xlabel('x position (meters)')

legend('1D','FDM - Coarse','FDM - Fine')    
    
    
%% Plotting qf,Ttip and nf along Lf/t
% 1D qf for plot along L/t
Tinf = 20; %Celsius
Tb = 85;%Celsius

h = 100; %W/m^2 K
k = 180; %W/m K

Lf = .012; % meter
L = Lf;
Lb = .003; % meter
W = 0.02;%in meters

N = 1;

t = Lf./[3  5  7  9 10];

P = 2.*(W+t);
Ac = t*W;
Af = 2.*((Lf.*W)+(Lf.*t));
At = W.^2 - N.*(W.*t) + N.*Af;

m = sqrt((h.*P)./(k.*Ac));

theta_b = Tb - Tinf;
M = sqrt(h.*P.*k.*Ac.*theta_b);
qf = M.*tanh(m.*L);

figure(8)
hold on

plot(Lf./t,qf)  

% FDM coarse qf for plot along L/t
fdmCqf = xlsread('fdmc.xlsx','B18:F18');
plot(Lf./t,fdmCqf) 

% FDM Fine qf for plot along L/t
fdmFqf = xlsread('fdmf.xlsx','B16:F16');
plot(Lf./t,fdmFqf) 


title('Fin Heat Transfer')
ylabel('qf (W/m^2)')
xlabel('Lf/t')

legend('1D','FDM - Coarse','FDM - Fine') 
hold off



%Ttip------------------------------------------------------
figure(9)

hold on

% FDM coarse qf for plot along L/t
fdmCTtip = xlsread('fdmc.xlsx','B23:F23');
plot(Lf./t,fdmCTtip) 

% FDM Fine qf for plot along L/t
fdmFqTtip = xlsread('fdmf.xlsx','B21:F21');
plot(Lf./t,fdmFqTtip) 



%nf------------------------------------------------------
figure(10)
