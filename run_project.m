clear all
close all

disp('Running model...')

model = l2('model');
%model.simulate(5, 'COM3');
model.simulate(80);
%model.plot();
                    %do not plot beliefs etc
model.plot({'anxiety','sitfac', 'hr',...
    'breathing_f', 'physical_state',...
    'chest_c','graph_bel_chest_pos',...
    'graph_bel_chest_trans','graph_bel_breathing_f',...
    'graph_breathing_f_error'});
%plot belief...


%model.export('');