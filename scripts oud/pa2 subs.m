%% deriving a formula from pairs of points
clc, clear all
syms a b c x1 x2 x3 y1 y2 y3 % a2 z
% f(x)...   f: y = ...
f1 = a*x1^2 + b*x1 + c % = y1
f2 = a*x2^2 + b*x2 + c % = y2
f3 = a*x3^2 + b*x3 + c % = y3

% Transcribe f1 and f2 to definitions of c (c1 and c2).
c1 = solve(f1 == y1, c) % c = y1 - a*x1^2 - b*x1
c2 = solve(f2 == y2, c)
% Equate c1 and c2 and derive a definition of a that is independent of c.
a1 = solve(c1 == c2, a)
% Substitute the original a in f2 and f3 with the new definition of a.
%This creates new definitions for f2 and f3 that are independent of the parameter a.
f3_no_a = subs(f3,a,a1)
f2_no_a = subs(f2,a,a1)

% Transcribe the new definitions f2 and f3 to definitions of b (b1 and b2).
b1 = solve(f3_no_a == y3, b)
b2 = solve(f2_no_a == y2, b)
% Equate b1 and b2 and derive a definition of c that is independent of
% both a and b. It consists of nothing more than the references to the values of x1-y3.
c_no_ab = solve(b1 == b2, c) % this c exist only of numbers

% Substitute the original c in f1 and the original f2 with the new 
% definition of c. This creates new definitions for f1 and f2 that are 
% independent of the parameter c.
f1 = subs(f1,c,c_no_ab)
f2 = subs(f2,c,c_no_ab)

% Transcribe the new definitions f1 and f2 to definitions of a (a1 and a2).
a1 = solve(f1 == y1, a)
a2 = solve(f2 == y2, a)
% Equate a1 and a2 and derive a definition of b that is independent of 
% both a and c. It consists of nothing more than the references to the 
%values of x1-y3.
b_no_ac = solve(a1 == a2, b)

% Substitute the original b in f1 with the new definition of b. This 
% creates new definitions for f1 that is independent of the parameters 
% a and c.
f1 = subs(f1,b,b_no_ac)
% Transcribe the new definitions f1 to a definitions of a. Now all the 
% definitions of each parameters is independent of the values of the other 
% parameters. 
a_no_bc = solve(f1 == y1, a)

a = simplify(a_no_bc)
% a=-(x1*y2 - x2*y1 - x1*y3 + x3*y1 + x2*y3 - x3*y2)/((x1 - x2)*(x1 - x3)*(x2 - x3))
b = simplify(b_no_ac)
% b=(x1^2*y2 - x2^2*y1 - x1^2*y3 + x3^2*y1 + x2^2*y3 - x3^2*y2)/((x1 - x2)*(x1 - x3)*(x2 - x3))
c = simplify(c_no_ab)
% c=-(- y3*x1^2*x2 + y2*x1^2*x3 + y3*x1*x2^2 - y2*x1*x3^2 - y1*x2^2*x3 + y1*x2*x3^2)/((x1 - x2)*(x1 - x3)*(x2 - x3))

% Example values for {x,y}
x1 =  90; y1 = 0.3;
x2 = 110; y2 = 0.45;
x3 = 125; y3 = 0.5;

% x1 = 70.1; y1 = 0.8230;
% x2 = 70.3; y2 = 0.8230;
% x3 = 70.2; y3 = 0.8230;
a =-(x1*y2 - x2*y1 - x1*y3 + x3*y1 + x2*y3 - x3*y2)/((x1 - x2)*(x1 - x3)*(x2 - x3))
b =(x1^2*y2 - x2^2*y1 - x1^2*y3 + x3^2*y1 + x2^2*y3 - x3^2*y2)/((x1 - x2)*(x1 - x3)*(x2 - x3))
c =-(- y3*x1^2*x2 + y2*x1^2*x3 + y3*x1*x2^2 - y2*x1*x3^2 - y1*x2^2*x3 + y1*x2*x3^2)/((x1 - x2)*(x1 - x3)*(x2 - x3))

% f(x) = a*x^2 + b*x + c
x=95;
y = a*x^2 + b*x + c
x=110;
y = a*x^2 + b*x + c
% dit werkt ook echt gewoon

