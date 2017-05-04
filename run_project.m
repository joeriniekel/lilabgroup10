clear all
close all

disp('Running model...')

model = l2('model');
%model.simulate(5, 'COM3');
model.simulate(50);
%model.plot();
                    %do not plot beliefs etc
model.plot({'anxiety','sitfac', 'hr', 'breathing_f', ...
  'physical_state','chest_c','graph_belief_chest_pos',...
  'belief'});
