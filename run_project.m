clear all
close all

disp('Running model...')

model = l2('lilab10');
%model.simulate(5, 'COM3');
model.simulate(50);
model.plot();
