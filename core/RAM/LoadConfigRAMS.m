function RamsData = LoadConfigRAMS(pFolder)

%% Loading main parameters

parSet = squeeze(struct2cell(ParseParams([pFolder 'MainRAMS.txt'])));


RamsData.freq = str2double(parSet(2,strcmp(parSet(1,:),'freq')));
RamsData.zs = str2double(parSet(2,strcmp(parSet(1,:),'zs')));
RamsData.zr = str2double(parSet(2,strcmp(parSet(1,:),'zr')));
RamsData.rmax = str2double(parSet(2,strcmp(parSet(1,:),'rmax')));
RamsData.dr = str2double(parSet(2,strcmp(parSet(1,:),'dr')));
RamsData.ndr = fix(str2double(parSet(2,strcmp(parSet(1,:),'ndr'))));
RamsData.zmax = str2double(parSet(2,strcmp(parSet(1,:),'zmax')));
RamsData.dz = str2double(parSet(2,strcmp(parSet(1,:),'dz')));
RamsData.ndz = fix(str2double(parSet(2,strcmp(parSet(1,:),'ndz'))));
RamsData.zmplt = str2double(parSet(2,strcmp(parSet(1,:),'zmplt')));
RamsData.c0 = str2double(parSet(2,strcmp(parSet(1,:),'c0')));
RamsData.np = fix(str2double(parSet(2,strcmp(parSet(1,:),'np'))));
RamsData.drProf = fix(str2double(parSet(2,strcmp(parSet(1,:),'drProf'))));
RamsData.dzHydr = 2;
RamsData.aModel = parSet(2,strcmp(parSet(1,:),'aModel'));


%% Various media parameters: bottom, hydrology, bathymetry

% bathymetry

RamsData.bath = dlmread([pFolder 'bath.txt']);

% geoacoustic parameters of the bottom

bData = importdata([pFolder 'bottom.txt']);
RamsData.bParams = bData.data(:,2:6);

% hydrology & bottom layers 

% if the file layers.txt exists, then profiles are constructed up to rmax
% if not, profiles are constructed up to the last hydrology profile
% the number of profiles corresponds to the cw size
% if bLayers has 1 column, then it is used for all profiles

if exist([pFolder 'layers.txt'],'file')
    RamsData.cw = HydroLoadRAMS([pFolder 'hydrology/'],RamsData.dzHydr, RamsData.drProf,RamsData.rmax);
    RamsData.bLayers = LayersLoadRAMS(pFolder, RamsData.drProf,RamsData.rmax);    
else
    RamsData.cw = HydroLoadRAMS([pFolder 'hydrology/'],RamsData.dzHydr, RamsData.drProf);
    RamsData.bLayers = bData.data(:,1);
end;

