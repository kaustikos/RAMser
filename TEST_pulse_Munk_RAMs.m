close all;
clear variables;
clc;
pFolder = 'TEST_Munk/';

% Load config and run the computation


% source function

dt = 10^(-4);
Tmax = 1;
tt = 0:dt:Tmax;

alph = 40;
bet = 100;
m = alph/bet;

%aa = (bet^2+alph^2)/alph;
phi = pi/2;
aa = (bet^2)*(m^2+1)/( m*cos(phi) - sin(phi)*m^2 );


ft_s = aa*exp(-bet*tt).*(-sin(alph*tt + phi) + m*cos(alph*tt + phi) + sin(phi) )  ;


figure;
plot(tt,ft_s,'linewidth',1.5);
grid on;
legend('dg(t)/dt');
xlabel('t, s');
xlim([0 0.1]);

dlmwrite([pFolder 'ft.txt'],[tt.' ft_s.'],'delimiter','\t','precision',8);

z = 0:2:4000;
zbar = 2*(z-1300)/1300;
eps = 0.00737;
c = 1500*(1 + eps*(zbar - 1 + exp(-zbar)));

dlmwrite([pFolder 'hydrology/0.hydr'],[z.' c.'],'delimiter','\t','precision',8);

figure;
plot(z,c,'linewidth',2);
grid on;       





zr = 2000;
R = 10000*(1:5);

%R = [50000 75000];
%R = [5000];
cmax = 1560;

% its OK! But use another value for c0 in RAMs

[t_r pt_r] = RunModelPulse(pFolder, zr, R,cmax);


for ii = 1:length(R)
    
    dlmwrite([pFolder int2str(R(ii)) '_3D_rams.txt'],[t_r(:,ii) pt_r(:,ii)],'delimiter','\t','precision',8);
       

    
end;


