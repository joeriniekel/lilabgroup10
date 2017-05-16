clc
disp('Running...')
clear all
close all
model = l2('model_arduino');tic;
n = 100;
%/dev/cu.usbmodem411 (Arduino/Genuine Uno)
%tty.usbmodem411
model.simulate(n,'COM1');
%a = arduino('/dev/tty.usbmodem411','Uno')
%model.simulate(5, 'COM3');
toc;
model.plot();

