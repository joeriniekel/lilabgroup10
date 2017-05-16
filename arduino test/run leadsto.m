disp('Running the leadsto tutorial model...')

clear all
close all
model = l2('01_leadsto_tutorial')
model.simulate(5, 'COM3');
model.plot();
