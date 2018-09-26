close all;
clear variables;
clc;
pFolder = 'wedge/';

% Load config and run the computation

RunModel(pFolder);

% plot resulting TL 

PlotResults(pFolder);