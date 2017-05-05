clear all
close all
clc
disp('Running model...')
model = l2('model');
%model.simulate(5, 'COM3');
model.simulate(50);
disp('Simulation finished');
%model.plot();
                    %do not plot beliefs etc
model.plot({'anxiety','sitfac', 'hr',...
    'breathing_f', 'physical_state',...
    'chest_c','graph_bel_chest_pos',...
    'graph_bel_breathing_f','graph_breathing_f_error'});
%plot belief...

%graph_breathing_f_error should generate 2 lines

%model.plot({'graph_bel_chest_trans'});
%legend('Undefined','Start of breathing cycle')

%model.plot({'graph_breathing_f_error'});
%legend('Deviation','Zero')
%model.export('../../plots/model.pdf',{'graph_breathing_f_error'});



