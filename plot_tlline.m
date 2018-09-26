set(0, 'DefaultAxesFontSize', 16, 'DefaultAxesFontName', 'Arial');
set(0, 'DefaultTextFontSize', 16, 'DefaultTextFontName', 'Arial'); 


figure;


% TL(r) plot at z = zr

% CoupleTL = dlmread([pFolder 'CoupleAtt.TL']); 
% 

hold all;

tlDat = dlmread([pFolder 'results/' 'tl.nLine_f360_zs40_zr100_100km.Txt']);
plot(tlDat(:,1)/1000,tlDat(:,3),'linewidth',2);

tlDat = dlmread([pFolder 'results/' 'tl.nLine_f360_zs35_zr100_100km.Txt']);
plot(tlDat(:,1)/1000,tlDat(:,3),'linewidth',2);

tlDat = dlmread([pFolder 'results/' 'tl.nLine_f360_zs30_zr100_100km.Txt']);
plot(tlDat(:,1)/1000,tlDat(:,3),'linewidth',2);

tlDat = dlmread([pFolder 'results/' 'tl.nLine_f360_zs25_zr100_100km.Txt']);
plot(tlDat(:,1)/1000,tlDat(:,3),'linewidth',2);

tlDat = dlmread([pFolder 'results/' 'tl.nLine_f360_zs20_zr100_100km.Txt']);
plot(tlDat(:,1)/1000,tlDat(:,3),'linewidth',2);

legend('40 m','35 m','30 m','25 m','20 m');
xlim([25 60]);
ylim([-160 -80]);
xlabel('r, km');
ylabel('TL, dB');

grid on;