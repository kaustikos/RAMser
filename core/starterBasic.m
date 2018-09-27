clc;
clear all;
model = BasicModel;
runner = BasicModelRunner;
launcher = Launcher(runner, model);
launcher.run();