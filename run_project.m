clear all
close all
clc
disp('Running model...')
model = l2('model');tic
n = 400;
global TRAINING TRAINING_BF TRAINING_HR
c1_bf = csvread('data/17-05 conditie 1 rust/bb_v69.csv');%499
c1_hr = csvread('data/17-05 conditie 1 rust/hr_v69.csv');
c2_bf = csvread('data/17-05 conditie 2 sport/bb_v81.csv');
c2_hr = csvread('data/17-05 conditie 2 sport/hr_v81.csv');
c3_bf = csvread('data/17-05 conditie 3 angst/bb_v32.csv');
c3_hr = csvread('data/17-05 conditie 3 angst/hr_v32.csv');
%als deze worden gebruikt worden de domein-waardes genegeerd
% dt = 0.18;
TRAINING = true;
TRAINING_BF = c1_bf;
TRAINING_HR = c1_hr;

% global RT_CHEST YY
% x = linspace(0,8);
% YY = sin(x);
% RT_CHEST = plot(x,YY);
% h.XDataSource = 'x';
% h.YDataSource = 'YY';
% linkdata on
global CHEST_Y1 RT_CHEST1 %CHEST_Y2 RT_CHEST2
CHEST_Y1 = zeros([1,100]); % CHEST_Y2 = [0]; subplot(1,2,2); 
RT_CHEST1 = stem(CHEST_Y1);linkdata on
% RT_CHEST1 = plot(CHEST_Y1);linkdata on
% subplot(1,2,1); RT_CHEST2 = plot(CHEST_Y2);linkdata on



%model.simulate(n, 'COM3');
model.simulate(n,'default','default');
disp('Simulation finished');
time = toc; disp('time/n');   disp(time/n);
%model.plot();
model.plot({...
    'chest_c'...
    ,'graph_bel_starting_dir'...
    ,'graph_bel_hr'...
    ,'graph_bel_breathing_f'...
    ,'graph_bel_anxiety','assessment'...
    ,'graph_bel_prev_ps','graph_original_hr'...
    'graph_breathing_f_diff'...
    ,'cycle_time'...
    ,'graph_des_starting_dir','support'
});
    %,'graph_breathing_f_error'...
    %'sitfac','anxiety_regulation','anxiety',...
    %'ps','hr','breathing_f','relative_c'...
    %,'starting_dir','performance'...
	%,'phase_shift'...
    
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



