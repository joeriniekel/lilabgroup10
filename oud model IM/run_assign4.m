clear all
close all
model = l2('assign4-L2')
model.simulate(28,'scenario2')
model.plot()
%model.plot({'feeling_of_loneliness','propose'})