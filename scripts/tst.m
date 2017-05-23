%%
prev_chest_c = 79.8;
avg_chest_c = 70;
range = 19;
relative_c = 0.4956;

phi = 0.77;
phi = 0.2475*pi;
f = 0.21;
dt = 0.5;

%plotting the breathing f
t = 0:0.001:10;
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


%plots:
t = 0:0.001:1;
%sinus curve
plot(t,sin(t*0.5*pi))
%inverse sinus
plot(t,asin(t)/pi * 2)


%%
clc
clear all
x = linspace(-100,100);
a = 0.001; b = 0.0031; c = 10;
f1 = a*x.^2 + b*x + c;
plot(x,f1)
syms x a b c
f1 = a*x.^2 + b*x + c
f2 = diff(f1)
% b + 2*a*x
f3 = b + 2*a*x;
%%
clear all
% x = linspace(-100,100);
% a = 0.001; b = 0.0031; c = 10;
% f3 = b + 2*a*x;
% find(min(f3) == f3)
% plot(x,f3)


mini = 30;
maxi = 250;
x = 30:2:250; %linspace
a = 0.001; b = 0.0031; c = 10;
f3 = b + 2*a*x;
min([1 1]);
min(f3)

% index = find(min(f2) == f2)
% f1(index)
% x(index)
% f2 = simplify(f2)
% pretty(f2)
% crit_pts = solve(f1)
%  -(b + (b^2 - 4*a*c)^(1/2))/(2*a)
%  -(b - (b^2 - 4*a*c)^(1/2))/(2*a)
%%
% function result = red_green_graph(model, trace, parameters, t)
   starting_dir = l2.getall(trace, t+1, 'desire', predicate('starting_dir', NaN)).arg{1}.arg{1};
   image=zeros(300,400,3);
   if starting_dir == '1 in'
	image(:,:,2)=1; %green
   elseif starting_dir == '3 out'
   	image(:,:,1)=1; %red
   else image(:,:,3)=1; %blue
   end
   figure, imshow(image)
% end

%%
image=zeros(300,400,3);
image(:,:,2)=1; %green
image(:,:,1)=1; %red
image(:,:,3)=1; %blue
figure, imshow(image)
%%
Y = [1, 5, 3;
     3, 2, 7;
     1, 5, 3;
     2, 6, 1];
h = area(Y,'LineStyle',':');
Y = [1 1];
h = area(Y,'LineStyle',':');
h(1).FaceColor = [1 0 0];% red
h(1).FaceColor = [0 0 1];% blue
h(1).FaceColor = [1 1 0];% yellow
h(1).FaceColor = [0 0.5 0];% green
% h(1).FaceColor = [0 0 0];% black
