clc;
clear all;
addpath('core')
addpath('tools')
addpath('RAM')
model = BasicModel;
runner = BasicModelRunner;
launcher = Launcher(runner, model);
launcher.run();