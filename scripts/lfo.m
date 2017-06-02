%%
close all, clear all, clc
p = 90;
q = 10;

x = linspace(0,2*pi * 3);
subplot(2,1,1); plot(x,p + q*sin(x))
hold on
subplot(2,1,1); plot(x,p + q*0.4*sin(0.7*x + 1))
legend('hr','modulator')
xlabel('time'); ylabel('hr');

subplot(2,1,2); plot(x,p + q*sin(0.6*x + 1) + q*sin(x))
legend('hr + modulator')
xlabel('time'); ylabel('hr');

