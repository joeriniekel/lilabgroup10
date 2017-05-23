clear all
close all
clc

% x = 1:20;
% y = rand(20,3);
% area(x,y)
% linkdata on

% y(10,:) = 0;

% p = 1:100
% q = p(1:10)
% plot(q)
% linkdata on

x = linspace(0,8);
y = sin(x);
% figure
h = plot(x,y);

h.XDataSource = 'x';
h.YDataSource = 'y';

y = sin(x.^3);
refreshdata

