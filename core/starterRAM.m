clc;
clear all;
model = BasicModel;
runner = RAMModelRunner;
launcher = Launcher(runner, model);
launcher.run();