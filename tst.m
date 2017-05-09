prev_chest_c = 79.8;
avg_chest_c = 70;
range = 19;
relative_c = 0.4956;

phi = 0.77;
phi = 0.2475*pi;
f = 0.21;
dt = 0.5;

%plotting the breathing f
t = 0:0.01:10;
f = 1.3;
plot(t,sin(t*f));grid on;
%the value of x(9) = sin(9*f)
t_9_f_9 = sin(9*f)
%assume that on t=9 the there is a change in f
f = 1.5;
plot(t,sin(t*f));grid on;
t_10_f_1_2 = sin(10*f)
%now the value of x seems to have skipped several time steps



%sin(2*pi* f * dt + phi)