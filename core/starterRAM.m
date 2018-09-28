clc;
clear all;
addpath('core')
addpath('tools')
addpath('RAM')

model = RAMModel;
model.rootFolder = 'RAM';
model.pFolder = [model.rootFolder, '/wedge/'];
runner = RAMModelRunner;
launcher = Launcher(runner, model);
launcher.run();
PlotResults(model.pFolder);