%%
clc; clear all; close all;

% length(c1_bf)
c1_bf = csvread('data/calculated/bb_c1_v61.csv');
c1_hr = csvread('data/calculated/hr_c1_v61.csv');
c2_bf = csvread('data/calculated/bb_c2_v51.csv');
c2_hr = csvread('data/calculated/hr_c2_v51.csv');



length(c2_bf)

from = 200;to = 310;
subplot(2,1,1);plot(c1_bf(from:to))
subplot(2,1,2);plot(c1_hr(from:to))
mean(c1_bf(from:to)) % 0.5159
mean(c1_hr(from:to)) % 67.2703

from = 250;to = 490;
subplot(2,1,1);plot(c2_bf(from:to))
subplot(2,1,2);plot(c2_hr(from:to))
mean(c2_bf(from:to)) % 0.6238
mean(c2_hr(from:to)) % 99.8270

%from = 60;to = 100;
%subplot(2,1,1);plot(c2_bf(from:to))
%subplot(2,1,2);plot(c2_hr(from:to))
%mean(c2_bf(from:to)) % 0.6237
%mean(c2_hr(from:to)) % 67.3010

from = 130;to = 180;
subplot(2,1,1);plot(c2_bf(from:to))
subplot(2,1,2);plot(c2_hr(from:to))
mean(c2_bf(from:to)) % 0.5224
mean(c2_hr(from:to)) % 75.1351
