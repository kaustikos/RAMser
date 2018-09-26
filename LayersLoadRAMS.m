function bLayers = LayersLoadRAMS(pFolder,drProf,rmax)

tmp = dlmread([pFolder 'layers.txt']);
rProfIn = tmp(1,:);


bLayersIn = tmp(2:end,:);
nLayers = size(bLayersIn,1);


nProf = fix(rmax/drProf)+1;



bLayers(1:nLayers,1:nProf) = interp1(rProfIn.', bLayersIn.' , drProf*(0:nProf-1).','linear','extrap').';    