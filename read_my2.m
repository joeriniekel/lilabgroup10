%%
clc; clear all; close all;

% length(c1_bf)
c1_bf = csvread('data/calculated/bb_c1_v61.csv');
c1_hr = csvread('data/calculated/hr_c1_v61.csv');
c2_bf = csvread('data/calculated/bb_c2_v51.csv');
c2_hr = csvread('data/calculated/hr_c2_v51.csv');
length(c2_bf)

%%
from = 1;to = 497;
subplot(2,1,1);plot(c1_bf(from:to))
subplot(2,1,2);plot(c1_hr(from:to))
%%
from = 50;to = 125;
subplot(2,1,1);plot(c1_bf(from:to))
subplot(2,1,2);plot(c1_hr(from:to))
mean(c1_bf(from:to)) % 0.4986
mean(c1_hr(from:to)) % 61.0997
%%
from = 260;to = 340;
subplot(2,1,1);plot(c1_bf(from:to))
subplot(2,1,2);plot(c1_hr(from:to))
mean(c1_bf(from:to)) % 0.4751
mean(c1_hr(from:to)) % 72.1079
%%
from = 55;to = 497;
subplot(2,1,1);plot(c2_bf(from:to))
subplot(2,1,2);plot(c2_hr(from:to))
%%
from = 270;to = 340;
subplot(2,1,1);plot(c2_bf(from:to))
subplot(2,1,2);plot(c2_hr(from:to))
mean(c2_bf(from:to)) % 0.5678
mean(c2_hr(from:to)) % 97.4585
%%
%from = 60;to = 100;
%subplot(2,1,1);plot(c2_bf(from:to))
%subplot(2,1,2);plot(c2_hr(from:to))
%mean(c2_bf(from:to)) % 0.6237
%mean(c2_hr(from:to)) % 67.3010
%%
from = 130;to = 180;
subplot(2,1,1);plot(c2_bf(from:to))
subplot(2,1,2);plot(c2_hr(from:to))
mean(c2_bf(from:to)) % 0.5224
mean(c2_hr(from:to)) % 75.1351
%%

from = 1;to = 497;
subplot(2,1,1);plot(c2_bf(from:to))
subplot(2,1,2);plot(c2_hr(from:to))
mean(c2_bf(from:to)) % 0.5224
mean(c2_hr(from:to)) % 75.1351
