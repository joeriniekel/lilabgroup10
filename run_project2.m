clear all; close all; clc
global N PLOT_BF HR_AXIS BF_AXIS
N = 498;
HR_AXIS = 40:200;   BF_AXIS = 0.01*HR_AXIS;

subplot(1,2,1);     PLOT_BF = plot(HR_AXIS,BF_AXIS);    axis([0 200 0 4]);
%xlabel='Heart Rate';ylabel='Breathing Frequency';Title='Relation between hr & bf';
linkdata on

global TRAINING TRAINING_BF TRAINING_HR
c1_bf = csvread('data/17-05 conditie 1 rust/bb_v69.csv');%499
c1_hr = csvread('data/17-05 conditie 1 rust/hr_v69.csv');
c2_bf = csvread('data/17-05 conditie 2 sport/bb_v81.csv');
c2_hr = csvread('data/17-05 conditie 2 sport/hr_v81.csv');
c3_bf = csvread('data/17-05 conditie 3 angst/bb_v32.csv');
c3_hr = csvread('data/17-05 conditie 3 angst/hr_v32.csv');
    % als deze worden gebruikt worden de domein-waardes genegeerd
    % dt = 0.18;
TRAINING    = false;
TRAINING_BF = c3_bf;    TRAINING_HR = c3_hr;

% global RT_CHEST YY
% x = linspace(0,8); YY = sin(x);
% RT_CHEST = plot(x,YY);
% h.XDataSource = 'x';  h.YDataSource = 'YY'; linkdata on
global CHEST_Y1 PLOT_CHEST1  %CHEST_Y2 RT_CHEST2
CHEST_Y1 = zeros([1,100]); % CHEST_Y2 = [0]; subplot(1,2,2); 
subplot(1,2,2);
PLOT_CHEST1 = stem(CHEST_Y1);linkdata on
% RT_CHEST1 = plot(CHEST_Y1);linkdata on
% subplot(1,2,1); RT_CHEST2 = plot(CHEST_Y2);linkdata on
% refreshdata(RT_CHEST1);


disp('Running model...')
model = l2('model');tic
%model = l2('model','COM3');tic
model.simulate(N, 'COM5');
%model.simulate(N,'default','default');
disp('Simulation finished');
time = toc; disp('time/N');   disp(time/N);
%model.plot();
model.plot({...
    'hr','breathing_f',...
    'chest_c',...
    'graph_bel_starting_dir'...
    'graph_bel_hr'...
    'graph_bel_breathing_f'...
    'graph_des_breathing_f'...
    'graph_breathing_f_diff'...
    'graph_bel_used_chest_range'...
    'graph_stable_hr','graph_stable_bf'...
    'graph_bel_anxiety','assessment'...
    'graph_bel_prev_ps','graph_original_hr'...
    'cycle_time'...
    'graph_des_starting_dir','support'...
    'adaption_2'
});
    %,'graph_breathing_f_error'...
    %'sitfac','anxiety_regulation','anxiety',...
    %'ps','hr','breathing_f','relative_c'...
    %,'starting_dir','performance'...
	%,'phase_shift'...
    
 %  ,'physical_state'...
 %   ,'graph_bel_breathing_acc','graph_bel_breathing_pattern'...

 

 
%graph_bel_breathing_pattern
%plot belief...

%graph_breathing_f_error should generate 2 lines

%model.plot({'graph_bel_chest_trans'});
%legend('Undefined','Start of breathing cycle')

%model.plot({'graph_breathing_f_error'});
%legend('Deviation','Zero')
%model.export('../../plots/model.pdf',{'graph_breathing_f_error'});



