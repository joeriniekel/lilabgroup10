clc

% length(c1_bf)
c1_bf = csvread('data/calculated/bb_c1_v56.csv');
c1_hr = csvread('data/calculated/hr_c1_v56.csv');
c2_bf = csvread('data/calculated/bb_c2_v44.csv');
c2_hr = csvread('data/calculated/hr_c2_v44.csv');
c3_bf = csvread('data/calculated/bb_c3_v69 long.csv');
c3_hr = csvread('data/calculated/hr_c3_v69 long.csv');

length(c1_bf)

from = 1;to = 399;
subplot(2,1,1);plot(c3_bf(from:to))
subplot(2,1,2);plot(c3_hr(from:to))
