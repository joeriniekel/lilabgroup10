c1_bf = csvread('data/17-05 conditie 1 rust/bb_v69.csv');%499
c1_hr = csvread('data/17-05 conditie 1 rust/hr_v69.csv');
c2_bf = csvread('data/17-05 conditie 2 sport/bb_v81.csv');
c2_hr = csvread('data/17-05 conditie 2 sport/hr_v81.csv');
c3_bf = csvread('data/17-05 conditie 3 angst/bb_v32.csv');
c3_hr = csvread('data/17-05 conditie 3 angst/hr_v32.csv');

c3_bf = c3_bf(400:899)
c3_hr = c3_hr(400:899)
subplot(3,2,1);plot(c1_bf)
subplot(3,2,2);plot(c1_hr)
subplot(3,2,3);plot(c2_bf)
subplot(3,2,4);plot(c2_hr)
subplot(3,2,5);plot(c3_bf)
subplot(3,2,6);plot(c3_hr)

subplot(3,1,1);plot(c1_bf)
subplot(3,1,2);plot(c2_bf)
subplot(3,1,3);plot(c3_bf)

% avg hr
mean(c1_hr)
mean(c2_hr)
mean(c3_hr)

mean(c1_bf)/mean(c1_hr)
mean(c2_bf)/mean(c2_hr)
mean(c3_bf)/mean(c3_hr)
