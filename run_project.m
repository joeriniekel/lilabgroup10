clear all
close all
clc
disp('Running model...')
model = l2('model');
%model.simulate(5, 'COM3');
model.simulate(180,'default','default');
disp('Simulation finished');
%model.plot();

model.plot({'sitfac','anxiety_regulation','anxiety',...
    'ps','hr','breathing_f','prev_relative_c'...
	,'prev_phase_shift','chest_c','starting_dir'...
    ,'graph_bel_starting_dir'...
    ,'graph_bel_breathing_f','graph_d_bf','graph_breathing_f_error'...
    ,'graph_d_hr'...
    ,'graph_bel_anxiety'...
    ,'graph_bel_prev_ps','graph_original_hr'...
});    

 %, 'physical_state'...
 %   ,'graph_bel_breathing_acc','graph_bel_breathing_pattern'...
    
%graph_bel_breathing_pattern
%plot belief...

%graph_breathing_f_error should generate 2 lines

%model.plot({'graph_bel_chest_trans'});
%legend('Undefined','Start of breathing cycle')

%model.plot({'graph_breathing_f_error'});
%legend('Deviation','Zero')
%model.export('../../plots/model.pdf',{'graph_breathing_f_error'});



