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
    'graph_bel_breathing_f'});
%plot belief...

model.plot({'graph_bel_chest_trans'});
legend('Undefined','Start of breathing cycle')

model.plot({'graph_breathing_f_error'});
legend('Deviation','Zero')
model.export('model.png',{'graph_breathing_f_error'});



