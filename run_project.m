clear all
close all
clc
disp('Running model...')
model = l2('model');tic
n = 50;
% global ard2;
% ard2 = 4;
%model.simulate(n, 'COM3');
model.simulate(n,'default','default');
disp('Simulation finished');
time = toc; disp('time/n');   disp(time/n);
%model.plot();
model.plot({'sitfac','anxiety_regulation','anxiety',...
    'ps','hr','breathing_f','relative_c'...
	,'phase_shift','chest_c','starting_dir'...
    ,'graph_bel_starting_dir'...
    ,'performance'...
    ,'graph_bel_breathing_f','graph_breathing_f_error'...
    ,'graph_bel_anxiety','assessment'...
    ,'graph_bel_prev_ps','graph_original_hr'...
    'graph_breathing_f_diff'...
    ,'cycle_time'...
    ,'graph_des_starting_dir','support'
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



