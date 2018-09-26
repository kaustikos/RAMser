close all;
clear variables;
clc;
pFolder = 'wedge_Pasha/';

% Load config and run the computation
% pFolder, zr, R, cmax
zr = 30;
R = 4000;
cmax = 1500;

RunModelPulse(pFolder, zr, R, cmax);

% plot resulting TL 

%PlotResults(pFolder);