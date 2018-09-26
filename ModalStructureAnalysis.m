function [mAmpl, mgv_ptb, mgv_fd, grAngles, Epercentil, nmmax, EpercentilZ] = ModalStructureAnalysis(pFolder,rranges,nmmax)

set(0, 'DefaultAxesFontSize', 16, 'DefaultAxesFontName', 'Arial');
set(0, 'DefaultTextFontSize', 16, 'DefaultTextFontName', 'Arial');


tlDat = dlmread([pFolder 'results/' 'tl.nLine.Txt']);

bottom(:,1) = tlDat(:,1);
bottom(:,2) = tlDat(:,end);


[aFieldTL, dr, dz, aFieldP] = ReadRamsBinary([pFolder 'results/']);

nr = size(aFieldTL,2);
nz = size(aFieldTL,1);
nrr = length(rranges);


r = (0:nr-1)*dr;
z = (0:nz-1)*dz;

hr = interp1(bottom(:,1),bottom(:,2),r);
hrr = interp1(bottom(:,1),bottom(:,2),rranges);

aFieldPrr(1:nrr,1:nz) = interp1(r,aFieldP.', rranges );

figure;
hold all;
for ii = 1:nrr
    plot(abs(aFieldPrr(ii,1:nz)),z,'linewidth',1.5);
end;
grid on;
ylabel('z, m');
xlabel('|P(z)|, Pa');
set(gca,'YDir','reverse');


%% modal decomposition

RamsData = LoadConfigRAMS(pFolder);

% !to change: bathData(end,2) = last bath. point
% only 1 bottom layer included

cw = 1450;
rhow = 1;

MP.LayersData = [[0 cw cw rhow rhow 0 0]; [RamsData.bath(end,2) cw RamsData.bParams(2,1) rhow RamsData.bParams(2,3) 0 0] ];

% for this case we only take the first hydro file, i.e. the ssp profile is
% the same for all "r" up to 70 km.

hydrFolder = [pFolder 'hydrology/'];
fHydrList = '*.hydr';
fHydrList = GetFiles([hydrFolder,fHydrList],'','ASC');
MP.HydrologyData = dlmread([hydrFolder fHydrList(end).name]);

cz = MP.HydrologyData;

opts.nmod = -1;
opts.Ngr = 3;
opts.Tgr = 3;

freq = RamsData.freq;



mAmpl(1:nrr,1:nmmax) = 0;

grAngles(1:nrr,1:nmmax) = 0;

mgv_ptb(1:nrr,1:nmmax) = 0;
mgv_fd(1:nrr,1:nmmax) = 0;

Epercentil(1:nrr,1:nmmax) = 0;

EpercentilZ(1:nrr,1:nz) = 0;

Eint(1:nrr) = 0;

for ii = 1:nrr
    disp('Processing range r=');
    disp(rranges(ii));
    
    MP.LayersData(2,1) = hrr(ii);
    
    if hrr(ii) <= 400;
        % shallow-water regime
        opts.Hb = hrr(ii)*6;
        dz0 = 0.5;
    else
        % deep-water regime
        opts.Hb = RamsData.zmax;
        dz0 = 2;
    end;
    
    [wnum, wmode] = ac_modesr(dz0,MP,freq,opts);
    %disp(wnum.');
    
    nzm = size(wmode,1);
    zm = dz0*(0:nzm-1);
    %nmc = size(wmode,2);
    
    nmc = length(wnum);
    
    wmode = wmode(:,1:nmc);
    
    wmode_g = interp1(zm,wmode,z,'linear','extrap');
    
    % inverse density for the scalar product
    
    gamma(1:nz) = 1;
    izb = find(z>=hrr(ii),1,'first');
    gamma(izb:nz) = 1/RamsData.bParams(2,3);
    
    iczh = find(cz(:,1)>=hrr(ii),1,'first');
    cmin = min(cz(1:iczh,2));
    kmax = 2*pi*freq/cmin;
    
    grAngles(ii,1:nmc) = 180*acos(wnum/kmax)/pi;
    
    mgv_ptb(ii,1:nmc) = ModesGroupVelocities(zm,freq,wnum,wmode,MP);
    
    df = 0.5;
    [wnum1, ~] = ac_modesr(dz0,MP,freq+df,opts);
    
    mgv_fd(ii,1:nmc) = 2*pi*df./(  wnum1(1:nmc) - wnum(1:nmc)  );
    
    % computing amplitudes via gamma-weighted scalar product
    
    Eint(ii) = dz*abs((aFieldPrr(ii,1:nz).*gamma)*( aFieldPrr(ii,1:nz)' ));
    
    mAmpl(ii,1:nmc) = (aFieldPrr(ii,1:nz).*gamma)*wmode_g(:,1:nmc)*dz;
    
    Eperc = 0;
    for jj = 1:nmc
        Eperc = Eperc + abs(mAmpl(ii,jj))^2;
        Epercentil(ii,jj) = Eperc/Eint(ii);
    end;
    
    EpercentilZ(ii,1:nz) = dz*cumtrapz((abs(aFieldPrr(ii,1:nz)).^2).*gamma,2);
    EpercentilZ(ii,1:nz) = EpercentilZ(ii,1:nz)/EpercentilZ(ii,nz);
    
end;



% TL(r,z) computed from model's TL output




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
title('z_s = 40 m,  f = 300 Hz');

