clear all
close all
clc
disp('Running model...')
model = l2('model');
%model.simulate(5, 'COM3');
model.simulate(40,'default','default');
disp('Simulation finished');
%model.plot();

model.plot({'anxiety','sitfac', 'hr'...
    ,'breathing_f', 'physical_state','chest_pos_phi'...
    ,'chest_c2','chest_c','graph_bel_chest_pos'...
    ,'graph_bel_breathing_f','graph_breathing_f_error'...
 %   ,'graph_bel_breathing_acc','graph_bel_breathing_pattern'...
    });
%graph_bel_breathing_pattern
%plot belief...

%graph_breathing_f_error should generate 2 lines

%model.plot({'graph_bel_chest_trans'});
%legend('Undefined','Start of breathing cycle')

%model.plot({'graph_breathing_f_error'});
%legend('Deviation','Zero')
%model.export('../../plots/model.pdf',{'graph_breathing_f_error'});



