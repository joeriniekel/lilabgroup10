disp('Running the leadsto tutorial model...')

clear all
close all
tic
model = l2('model2');
global n;
n = 500;%max 5000
model.simulate(n, 'COM5'); %COM3 COM4 COM5
time = toc;
disp(time/n)
model.plot();

