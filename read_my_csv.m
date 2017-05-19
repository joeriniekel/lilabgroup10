clear all
close all
clc

% length(c1_bf)
c1_bf = csvread('data/17-05 conditie 1 rust/bb_v69.csv');
c1_hr = csvread('data/17-05 conditie 1 rust/hr_v69.csv');
c2_bf = csvread('data/17-05 conditie 2 sport/bb_v81.csv');
c2_hr = csvread('data/17-05 conditie 2 sport/hr_v81.csv');
c3_bf = csvread('data/17-05 conditie 3 angst/bb_v32.csv');
c3_hr = csvread('data/17-05 conditie 3 angst/hr_v32.csv');
length(c1_bf)

from = 1;to = 100;
subplot(2,1,1);plot(c1_bf(from:to))
subplot(2,1,2);plot(c1_hr(from:to))
