close all;
clear all;
clc;
pFolder = 'TEST_flat/';

set(0, 'DefaultAxesFontSize', 16, 'DefaultAxesFontName', 'Arial');
set(0, 'DefaultTextFontSize', 16, 'DefaultTextFontName', 'Arial');

% LOAD PULSE DATA
xx = dlmread([pFolder 'ft.txt']);
dt = xx(2,1)-xx(1,1);
ft_s = xx(:,2); 
% END PULSE DATA

tau = 0.01;
sig = 0.003;

t = 0:dt:0.5;
ft = (10^3)*(t-tau).^1.*exp( -(t-tau).^2/(sig^2) );
ft_s = ft;

ampl_tol = 0.1;


% LOAD RAMS
RamsData = LoadConfigRAMS(pFolder);
cmax = max( RamsData.bParams(:,1) );
% END LOAD RAMS

nf = length(ft_s);
Tf = dt*(nf-1);
fs = 1/dt;

figure;
plot( dt*(0:nf-1),ft_s,'linewidth',1.5 );
xlabel('t, s');
ylabel('p(t), s');
title('Source function (emitted signal)');
grid on;

NFFT = 2^nextpow2(nf);

Spf = ifft(ft_s,NFFT);
f = fs/2*linspace(0,1,NFFT/2+1);
ampl_sp = max(abs(Spf));
ind_f_last = find( abs( Spf(1:NFFT/2+1) ) > ampl_sp*ampl_tol,1,'last');
disp('Upper bound of freq spectrum');
disp(f(ind_f_last) );
ind_f_first = find( abs( Spf(1:NFFT/2+1) ) > ampl_sp*ampl_tol,1,'first');
disp('Lower bound of freq spectrum');
disp(f(ind_f_first) );

figure;
plot(f, abs(Spf(1:NFFT/2+1)),'linewidth',2,'color','black' );
grid on;
xlabel('f, Hz');
title('Frequency spectrum of the source function','linewidth',1.5);
xlim([0 f(ind_f_last)]);



dlmwrite('ft.txt',[t.' ft.'],'delimiter','\t','precision',8);
