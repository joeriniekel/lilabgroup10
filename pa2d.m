%% deriving a formula from pairs of points
clc, clear all
syms a b c x1 x2 y1 y2 y3 lhr
% f(x)...   f: y = ...
f1 = b*(x1-lhr) + c % = y1
f2 = b*(x2-lhr) + c % = y2
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

% done

x1 =  90; y1 = 0.9;
x2 = 110; y2 = 1;

x1 =  67.2703; y1 = 0.5159;
x2 = 99.8270; y2 = 0.6238;
%x3 = 75.1351; y3 = 0.5224;

lhr = 40;
%a = 4.2857e-04 % when filled in in the 'advanced' formula
%b = (y1 - y2 - a*(lhr - x1)^2 + a*(lhr - x2)^2)/(x1 - x2)
%c = y1 - a*(lhr - x1)^2 + ((lhr - x1)*(y1 - y2 - a*(lhr - x1)^2 + a*(lhr - x2)^2))/(x1 - x2)
b = (y1 - y2)/(x1 - x2)
c = y1 + ((lhr - x1)*(y1 - y2))/(x1 - x2)

%b =    0.0033
%c =    0.4255

x = linspace(0,200);
plot(x,b*(x-lhr) + c)





