close all;
clear variables;
clc;
pFolder = 'TEST_flat/';

% Load config and run the computation

zr = 50;
R = [5000 7000 8000];

cmax = 1700;
[t_r pt_r] = RunModelPulse(pFolder, zr, R, cmax);

for ii = 1:length(R)
    
    dlmwrite([pFolder int2str(R(ii)) '_rams.txt'],[t_r(:,ii) pt_r(:,ii)],'delimiter','\t','precision',8);
    
end;

% xx = dlmread([pFolder '8000_ac_modes.txt']);
% 
% xx = dlmread([pFolder 'ft.txt']);
% yy = dlmread([pFolder '8000_rams.txt']);
% 
% figure;
% hold all;
% plot(xx(:,1),xx(:,2));
% grid on;
% xlim([0 0.04]);
% xlabel('t, s');
% ylabel('p(t), Pa');


%plot(yy(:,1),yy(:,2));