% clear all
% close all
% clc
x = linspace(50,250)
% x2 = x.^2;
% p = plot(x2);linkdata on
% time = 0.2;
% pause(time);

% x2 = x*5;       refreshdata(p);pause(time);
% x2 = x.^0.5;    refreshdata(p);pause(time);
% x2 = x.^2 + 5;  refreshdata(p);pause(time);
% x2 = x.^2 + 5*x;refreshdata(p);pause(time);
% j=0;
% for i=0:0.1:5
%    x2 = j*x.^3 + i*x.^2;
%    j = j * -1 + 0.1;
%    refreshdata(p);pause(time);
% % end
% subplot(1,2,1); plot(x.^2 + x)
% % subplot(1,2,2); plot(2.^x + x)
% subplot(1,2,2); plot(x.^2 - 0.115*x.^3)
% plot(x,x.^2 - 0.115*x.^3)

% fig = figure(1);
% set(fig, 'Position', [10 10 10 10])
% figure('units','normalized','position',[.1 .1 .4 .4])

plot(x, 0.1*x - 5);
plot(x, 0.2*x.^2 - 20*x);

axis([50 250 0 100])
% axis([0 2*pi -1.5 1.5])

% bf = a*hr + b*hr^2....
bf = 0.3;
hr = 90;
a = bf/hr;
b = bf/(hr^2);
a
b

%%
syms x y
y = x*5 + 2;
double(subs(y,x,3))
%%
syms x y z a b
y = a * x
% subs({a,b},{x,y},{1,2})
subs(y,{x,a},{2,5})

%% matrices

syms a b c
M = [a b c; b c a; c a b]
M = [a b c; b c a; c a b]
M = subs(M, b, a + 1)
M = subs(M, M(1,3), a + 2)
[V,E] = eig(M)

%% trying
% stage 1; 1 var and 1 pair of coordinates
syms x1 y1 x2 y2 a b
% y1 = a*x1 + b
% b = y1 / (a*x1)
b = y1 / x1
ans = subs(b,{x1,y1},{2,5}) %double(ans)

%stage 2: 2 pairs
syms x1 y1 x2 y2 a b
% y1 = a*x1 + b
% y2 = a*x2 + b
% a = (y1 - b) / x1
% a = (y2 - b) / x2
% b = y1 / (a*x1)
% b = y2 / (a*x2)
% substitute a with a = (y1 - b) / x1
% b = y2 /(  ((y1 - b) / x1)   * x2)
% omschrijven:
% b = y2 /(  ((y1 - b) / x1)   * x2)

%poging 2
% a = (y1 - b) / x1   =   a = (y2 - b) / x2
% (y1 - b) / x1 = (y2 - b) / x2
% y1 - b = (y2 - b) / x2    * x1
% -b = (y2 - b) / x2    * x1  - y1
% b = - ((y2 - b) / x2    * x1  - y1)
%   %  - ((y2 - b) / x2 == -(y2/x2 - b/x2) = -y2/x2 + -b/x2
% b = (-y2/x2 + -b/x2)    * x1  - y1)
% b = -y2*x1/x2 + -b*x1/x2 - y1
% b - (-b*x1/x2)          = -y2*x1/x2 - y1
% b  + b*x1/x2            = -y2*x1/x2 - y1
% b * x2/x2  + b*x1/x2    = -y2*x1/x2 - y1
% b*x2/x2  + b*x1/x2  =   (b*x2 + b*x1) /x2   = -y2*x1/x2 - y1
% b*x2 + b*x1         = x2 * (-y2*x1/x2 - y1)
% b(x2 + x1)          = x2 * (-y2*x1/x2 - y1)
b = (x2 * (-y2*x1/x2 - y1))  /   (x2 + x1)

ans_b = subs(b,{y2,x2,x1,y1},{2,5,4,1})
ans_b = double(ans_b);

a = (y1 - b)
subs(a,{y1,b},{1,ans_b})


% y2 = a*x2 + b
% M = []
% subs(M, 
% [V,E] = eig(M)


%% test

syms x1 y1 x2 y2 a b z
y1 = a*x1 + b
y2 = a*x2 + b
% a = (y2 - b) / x2
z = (y2/x2 - b/x2) % z = ((y2 - b) / x2) does not work for some reason
subs(y1,a,z)

%% test2

syms a b x1 y1 x2 y2 z
% y1 = a*x1 + b
% y2 = a*x2 + b
% a = (y1-b)/x1
% a = (y2-b)/x2
ans = sym('(y1-b)/x1 = (y2-b)/x2')
% ans = sym('a + x1 == a / y1')
ans_to_b = solve(ans, b) % b = (x1*y2 - x2*y1)/(x1 - x2)
a = (y1-b)/x1
ans_to_a = subs(a,b,ans_to_b) % a = (y1 - (x1*y2 - x2*y1)/(x1 - x2))/x1
% if 2 pairs of {x,y} are known,
% a and b can be derived

%% Finally
clc
syms a b x1 y1 x2 y2 z
y1 = a*x1 + b
y2 = a*x2 + b
% b = y1 - a*x1
% b = y2 - a*x2
% b == b
comparison = sym('y1 - a*x1 = y2 - a*x2') % equation
ans_to_a = solve(comparison,a)
ans_to_b1 = subs(y1,a,ans_to_a)
% if 2 pairs of {x,y} are known,
% a and b can be derived

% % controle1:
% ans_to_b2 = subs(y2,a,ans_to_a)
% % controle2:
% % a = (y1-b)/x1 = (y2-b)/x2
% comparison2 = sym('(y1-b)/x1 = (y2-b)/x2')
% ans_to_b3 = solve(comparison,a)
% ans_to_a2 = subs(y1,b,ans_to_b3)
% % controle3:
% ans_to_a3 = subs(y2,b,ans_to_b3)


%% matrix test
syms a b x1 y1 x2 y2 z

M = [a*x1 b; a*x2 b]
% [y1 y2]
 
% M(2,2)

% [V,E] = eig(M)

%% next step: stage 3
clc
syms a b c x1 x2 x3 y1 y2 y3 a2 z

y1 = a*x1^2 + b*x1 + c
y2 = a*x2^2 + b*x2 + c
y3 = a*x3^2 + b*x3 + c

a1 = solve(y1,a) % y1 - (b*x1 + c*x1^2)
a2 = solve(y2,a) % y2 - (b*x2 + c*x2^2)
equation1 = sym(a1 == a2)
b_a1a2 = solve(equation1,b)
c_a1a2 = solve(equation1,c)

% y3_without_b = subs(y3,b,b_a1a2)
% a_without_b = solve(y3_without_b,a)
% c3 = solve(y3_without_b,c)





% a4 = subs(a1,c,c3)
% equation = sym(a == a4)
% a5 = solve(equation,a)
% 





% a1 = solve(y1,a) % y1 - (b*x1 + c*x1^2)
% a2 = solve(y2,a) % y2 - (b*x2 + c*x2^2)
% % a3 = solve(y3,a) % y3 - (b*x3 + c*x3^2)
% 
% equation1 = sym(a1 == a2)
% % equation2 = sym(a1 == a3)
% % equation3 = sym(a2 == a3)
% 
% b_a1a2 = solve(equation1,b)
% 
% y1b = subs(y1,b,b_a1a2)
% y2b = subs(y2,b,b_a1a2)
% y3b = subs(y3,b,b_a1a2)
% 
% a1 = solve(y1b,a)
% a2 = solve(y2b,a)
% a3 = solve(y3b,a)
% 
% equation2 = sym(a1 == a3)
% cxxx = solve(equation2,c)




% b_12 = solve(equation1,b)
% b_13 = solve(equation2,b)
% b_23 = solve(equation3,b)
% 
% c1 = solve(b_12,c)
% c1 = solve(b_13,c)
% c1 = solve(b_23,c)



% % % % 
% % % % y1_without_a = subs(y1,a,a2)
% % % % y2_without_a = subs(y2,a,a3)
% % % % %no y3
% % % % % y3_without_a = subs(y3,a,a1) %x2, x3
% % % % 
% % % % b1 = solve(y1_without_a,b)
% % % % b2 = solve(y2_without_a,b)
% % % % %
% % % % b3 = solve(y3,b)
% % % % 
% % % % % y1_without_ab = subs(y1_without_a,b,b1)
% % % % y2_without_ab = subs(y2_without_a,b,b3)
% % % % c1 = solve(y1_without_ab,c) % c= 0
% % % % c2 = solve(y2_without_ab,c)

% 
% equation = sym(b1 == b3)
% equation = sym(b1 == b3)
% a_alt = solve(equation,a)
% c_alt = solve(equation,c)



% y3_without_b = subs(y1,b,b3) %x2, x3
% 
% b1 = solve(y1_without_a,b)
% b2 = solve(y2_without_a,b)
% % no b3
% b3 = solve(y3_without_a,b)

% y1_without_ab = subs(y2_without_a,b,b1)
% y2_without_ab = subs(y2_without_a,b,b1)


% c1 = solve(y3_without_b,c)
% c2 = solve(y2_without_ab,c)

%% ...
syms x
f = x^3 - 15*x^2 - 24*x + 350
A = magic(3);
b = sym2poly(f)

syms a b x
f = 2*x^2 + 2*x
b = sym2poly(f)

%% xxxxxxx
clc
syms a b c x1 x2 x3 y1 y2 y3 a2 z
% f(x)...   f: y = ...
f1 = a*x1^2 + b*x1 + c % = y1
f2 = a*x2^2 + b*x2 + c % = y2
f3 = a*x3^2 + b*x3 + c % = y3

c1 = solve(f1 == y1, c) % c = y1 - a*x1^2 - b*x1
c2 = solve(f2 == y2, c)
a1 = solve(c1 == c2, a)

f3_no_a = subs(f3,a,a1)
f2_no_a = subs(f2,a,a1)

b1 = solve(f3_no_a == y3, b)
b2 = solve(f2_no_a == y2, b)
c_no_ab = solve(b1 == b2, c) % this c exist only of numbers

f1 = subs(f1,c,c_no_ab)
f2 = subs(f2,c,c_no_ab)

a1 = solve(f1 == y1, a)
a2 = solve(f2 == y2, a)
b_no_ac = solve(a1 == a2, b)

f1 = subs(f1,b,b_no_ac)
a_no_bc = solve(f1 == y1, a)

a = simplify(a_no_bc)
b = simplify(b_no_ac)
c = simplify(c_no_ab)


% a2 = solve(y2,a) % y2 - (b*x2 + c*x2^2)
% equation1 = sym(a1 == a2)
% b_a1a2 = solve(equation1,b)
% c_a1a2 = solve(equation1,c)



