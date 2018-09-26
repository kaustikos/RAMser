function [aFieldTL]=PlotResults(pFolder)

set(0, 'DefaultAxesFontSize', 16, 'DefaultAxesFontName', 'Arial');
set(0, 'DefaultTextFontSize', 16, 'DefaultTextFontName', 'Arial'); 


tlDat = dlmread([pFolder 'results/' 'tl.nLine.Txt']);

% TL(r) plot at z = zr

% CoupleTL = dlmread([pFolder 'CoupleAtt.TL']); 
% 

figure;
plot(tlDat(:,1)/1000,tlDat(:,3),'linewidth',2,'color','black');
% hold on;
% plot(CoupleTL(:,1),-CoupleTL(:,2),'color','red','linestyle','--','linewidth',2);


bottom(:,1) = tlDat(:,1);
bottom(:,2) = tlDat(:,end);


[aFieldTL, dr, dz, aFieldP] = ReadRamsBinary([pFolder 'results/']);

nr = size(aFieldTL,2);
nz = size(aFieldTL,1);

% TL(r,z) computed from model's TL output



figure;
imagesc((0:nr-1)*dr,(0:nz-1)*dz,aFieldTL);
caxis([-90 -40]);
hold on;
plot(bottom(:,1),bottom(:,2),'linewidth',2,'color','black')

% TL(r,z) computed from the complex pressure field
% IMPORTANT! scaling of the field by 4*pi is already applied!



figure;
imagesc((0:nr-1)*dr/1000,(0:nz-1)*dz,20*log10(abs(aFieldP)));
caxis([-130 -80]);
hold on;
plot(bottom(:,1)/1000,bottom(:,2),'linewidth',2,'color','white');
plot(bottom(:,1)/1000,bottom(:,2),'linewidth',1,'color','black');
xlabel('r, km');
ylabel('z, m');
colorbar;
%title('z_s = 35 m,  f = 490 Hz');

