function [wNum, mAmpl] = ModeDecomposition(pFolder,nmod,dz)

RamsData = LoadConfigRAMS(pFolder);


% !to change: bathData(end,2) = last bath. point
% only 1 bottom layer included

MP.LayersData = [[0 1450 1450 1 1 0 0]; [RamsData.bath(end,2) 1450 RamsData.bParams(1) 1 RamsData.bParams(3) 0 0] ];

% does it really work? does it take the last file?

hydrFolder = [pFolder 'hydrology/'];
fHydrList = '*.hydr';
fHydrList = GetFiles([hydrFolder,fHydrList],'','ASC');
MP.HydrologyData = dlmread([hydrFolder fHydrList(end).name]);

% loading RAMS computations results

[~, ~, dzf, aFieldP] = ReadRamsBinary([pFolder 'results/']);


% z grid in RAMS is usually to fine for spectral problem,
% we take a somewhat coarser one

nzf = size(aFieldP,1);
zf = (0:nzf-1)*dzf;

% TL(r,z) computed from the complex pressure field
% IMPORTANT! scaling of the field by 4*pi is already applied!


% the last slice P(rmax,z)

aFieldPzf = (aFieldP(1:nzf,end)).';

% new grid: magic number!

dz = 4;
z = 0:dz:zf(end);
nz = length(z);
aFieldPz = interp1(zf,aFieldPzf,z,'linear','extrap');

% inverse density for the scalar product

gamma(1:nz) = 1;
izb = find(z>=RamsData.bath(end,2),1,'first');
gamma(izb:nz) = 1/RamsData.bParams(3);

% solving spectral problem

% to use Richardson extrap. we have to apply correction to eigenfunctions!

[wNum, wmode] = ac_modesr(z,MP,RamsData.freq,nmod);

% computing amplitudes via gamma-weighted scalar product

mAmpl = (aFieldPz.*gamma)*wmode*dz;






