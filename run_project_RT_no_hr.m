clear all; close all; clc
global N REAL_TIME_INPUT TRAINING SOUND SAVE_DATA LIMIT_DT

% -------------------------
% CONFIG
% -------------------------

N               = 150;     % number of timesteps to render
SOUND           = true;    % use audio feedback for support
REAL_TIME_INPUT = true;   % use realtime input data
TRAINING        = false;   % use previously generated input from csv
SAVE_DATA       = true;   % save bf and hr to csv files
LIMIT_DT        = 0.6;

% -------------------------
% PLOTTING - realtime
% -------------------------

global PLOT_BF HR_AXIS BF_AXIS
HR_AXIS = 40:200;   BF_AXIS = 0.01*HR_AXIS;
subplot(2,2,1);     PLOT_BF = plot(HR_AXIS,BF_AXIS);    axis([0 200 0 4]);
%xlabel='Heart Rate';ylabel='Breathing Frequency';Title='Relation between hr & bf';
linkdata on


% global RT_CHEST YY
% x = linspace(0,8); YY = sin(x);
% RT_CHEST = plot(x,YY);
% h.XDataSource = 'x';  h.YDataSource = 'YY'; linkdata on
global CHEST_Y1 PLOT_CHEST1  %CHEST_Y2 RT_CHEST2
CHEST_Y1 = zeros([1,100]); % CHEST_Y2 = [0]; subplot(1,2,2); 
subplot(2,2,2);
PLOT_CHEST1 = stem(CHEST_Y1);linkdata on
% RT_CHEST1 = plot(CHEST_Y1);linkdata on
% subplot(1,2,1); RT_CHEST2 = plot(CHEST_Y2);linkdata on
% refreshdata(RT_CHEST1);

global PLOT_COLOR PLOT_TXT
subplot(2,2,3); PLOT_COLOR = area([1 1]);
PLOT_COLOR(1).FaceColor = [0 0 0];% red
% PLOT_COLOR(1).FaceColor = [0 0 1];% blue
subplot(2,2,4); PLOT_TXT = 'data/none.jpg'; imshow(PLOT_TXT); linkdata on
% imshow('data/none.jpg')



% -------------------------
% RUN
% -------------------------

disp('Running model...'); 
if REAL_TIME_INPUT && TRAINING, TRAINING = false; disp('-- WARNING not TRAINING --'); end;
if TRAINING, disp('using .csv data'); end;
model = l2('model_rt_no_hr');
tic;

% model.simulate(N, 'COM5');
model.simulate(N,'COM5','default','default');

disp('Simulation finished');
toc; % time = toc; disp('time/N');   disp(time/N);




% -------------------------
% PLOT full simulation
% -------------------------

model.plot({...
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



