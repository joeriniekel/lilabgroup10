%% sine and altered sine
close all, clear all, clc
x = linspace(0,1.2);
line([1 1],get(axes,'YLim'),'Color',[1 0 0]); hold on
f1 = 2.1;
phi = asin(sin(2*pi*1*f1));

f2 = 1.5;
plot(x,0.5 + 0.5*sin(2*pi*x*f1)); hold on
plot(x,0.5 + 0.5*sin(2*pi*(x+1)*f2 + phi)); hold on
axis([0 1.5 0 1])
legend('t2=0','original sine','altered sine')

%% lfo/modulator
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

