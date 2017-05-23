%%
clc; clear all; close all;

% length(c1_bf)
c1_bf = csvread('data/calculated/bb_c1_v14.csv');
c1_hr = csvread('data/calculated/hr_c1_v14.csv');
% c2_bf = csvread('data/calculated/bb_c2_v86.csv');
% c2_hr = csvread('data/calculated/hr_c2_v86.csv');
% c3_bf = csvread('data/calculated/bb_c3_v69 long.csv');
% c3_hr = csvread('data/calculated/hr_c3_v69 long.csv');



length(c1_bf)

from = 1;to = 497;
subplot(2,1,1);plot(c1_bf(from:to))
subplot(2,1,2);plot(c1_hr(from:to))

%mean(c2_bf(70:170)) % 0.7229
%mean(c2_hr(70:170)) % 69.7205

%mean(c2_bf(270:310)) % 0.9096
%mean(c2_hr(270:310)) % 97.9880



