%% deriving a formula from pairs of points
clc, clear all
syms a b c x1 x2 y1 y2 y3 lhr
% f(x)...   f: y = ...
f1 = a*(x1)^2 + b*(x1) + c % = y1
f2 = a*(x2)^2 + b*(x2) + c % = y2
% a is fixed

% Transcribe f1 and f2 to definitions of c (c1 and c2).
c1 = solve(f1 == y1, c) % c = y1 - a*x1^2 - b*x1
c2 = solve(f2 == y2, c)
% Equate c1 and c2 and derive a definition of b that is independent of c.
b_no_c = solve(c1 == c2, b)
% Substitute b with the new definition of b to derive an update formula 
% for f1.
f1_2 = subs(f1,b,b_no_c)

% Transcribe the new formula for f1 to a new definition of c (c3).
c3 = solve(f1_2 == y1, c)

b = b_no_c
c = c3

% done ------------

x1 =  90; y1 = 0.9;
x2 = 110; y2 = 1;
%x3 = 125; y3 = 1.3; % unused for now
x1 = 58.0938; y1 = 0.4995;
x2 = 71.29306667; y2 = 0.495;
x3 = 94.14; y3 = 0.6149;

% lhr = 40;
a = 1.5505e-04 % when filled in in the 'advanced' formula
b = (- a*x1^2 + a*x2^2 + y1 - y2)/(x1 - x2)
c = (a*x1^2*x2 - a*x1*x2^2 + y2*x1 - y1*x2)/(x1 - x2)


x = linspace(0,200);
plot(x,a*(x).^2 + b*(x) + c)
%%
a =   1.5505e-04
b =   -0.0204
c =    1.1615

x = linspace(0,200);
plot(x,a*(x).^2 + b*(x) + c)

