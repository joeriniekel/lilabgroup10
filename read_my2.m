%%
clc; clear all; close all;

% length(c1_bf)
c1_bf = csvread('data/calculated/bb_c1_v61.csv');
c1_hr = csvread('data/calculated/hr_c1_v61.csv');
c2_bf = csvread('data/calculated/bb_c2_v51.csv');
c2_hr = csvread('data/calculated/hr_c2_v51.csv');
length(c2_bf)


%% 1

%%
from = 1;to = 497;
subplot(2,1,1);plot(c1_bf(from:to))
subplot(2,1,2);plot(c1_hr(from:to))
%%
from = 50;to = 120;
subplot(2,1,1);plot(c1_bf(from:to))
subplot(2,1,2);plot(c1_hr(from:to))
mean(c1_bf(from:to)) % 0.4995
mean(c1_hr(from:to)) % 61.0938
%%
from = 260;to = 310;
subplot(2,1,1);plot(c1_bf(from:to))
subplot(2,1,2);plot(c1_hr(from:to))
mean(c1_bf(from:to)) % 0.5102
mean(c1_hr(from:to)) % 71.8111


%% 2

%%
from = 1;to = 497;
subplot(2,1,1);plot(c1_bf(from:to))
subplot(2,1,2);plot(c1_hr(from:to))
%%
from = 50;to = 120;
subplot(2,1,1);plot(c2_bf(from:to))
subplot(2,1,2);plot(c2_hr(from:to))
mean(c2_bf(from:to)) % 0.4786
mean(c2_hr(from:to)) % 67.9941
%%
from = 170;to = 180; % (kort)
subplot(2,1,1);plot(c2_bf(from:to))
subplot(2,1,2);plot(c2_hr(from:to))
mean(c2_bf(from:to)) % 0.4962
mean(c2_hr(from:to)) % 74.0740 %(kort)
%%
from = 220;to = 310;
subplot(2,1,1);plot(c2_bf(from:to))
subplot(2,1,2);plot(c2_hr(from:to))
mean(c2_bf(from:to)) % 0.6149
mean(c2_hr(from:to)) % 94.1400


