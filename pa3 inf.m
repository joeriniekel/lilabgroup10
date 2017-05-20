%% middag

clc
% bf = 
%     hr + a
%     a*hr + b 
%     a*hr^2 + b*hr + c
%     % etc

% % Thus
% % b = 1;
% % step 1
%     % bf = hr + a
% a = bf - b*hr
% %step 2
%     % bf = b*hr + a
% b = (bf-a)/hr
% %step 3
% verder met step 1
% bf = b*hr + a
% a = bf - b*hr

%formula will be bf = 0.0313 * hr -1.5536
b = 1;
a = 1;
hr_plot = 50:250;
bf_plot = a*hr_plot+b;
p = plot(hr_plot,bf_plot);linkdata on
time = 0.1;

for i=1:100
%t=1
hr = 100;   bf = 1.5764;
a = bf - b*hr;
bf_plot = a*hr_plot+b;
refreshdata(p)
pause(time)
%t=2
hr = 101;   bf = 1.6077;
b = (bf-a)/hr;
bf_plot = a*hr_plot+b;
refreshdata(p)
pause(time)
%t=3
hr = 90;   bf = 1.2634;
a = bf - b*hr;
bf_plot = a*hr_plot+b;
refreshdata(p)
pause(time)
%t=4
hr = 98;   bf = 1.5138;
b = (bf-a)/hr;
bf_plot = a*hr_plot+b;
refreshdata(p)
pause(time)
end
% repeat this many times and it will lead to the same result.

%% stage 3
clc
% bf = 
%     hr + a
%     a*hr + b 
%     a*hr^2 + b*hr + c
%     % etc

% % Thus
%formula will be bf = -1.1905e-04 * hr^2 + 0.0313*hr -1.5536 

b = 1;a = 1;c = 1;
hr_plot = 50:250;
bf_plot = a*hr_plot.^2+b*hr_plot+c;
p = plot(hr_plot,bf_plot);linkdata on
time = 0.1;

for i=1:100
    %t=1
    hr = 100;   bf = 0.3859;
        % bf = a*hr^2 + b*hr + c
        % (bf - b*hr - c) = a*hr^2
    a = (bf - b*hr - c) / hr;

    bf_plot = a*hr_plot.^2+b*hr_plot+c;
    refreshdata(p)
    pause(time)
    %t=2
    hr = 100;   bf = 0.4881;
        % bf = a*hr^2 + b*hr + c
    b = (bf - a*hr^2 - c)/hr;
    
	bf_plot = a*hr_plot.^2+b*hr_plot+c;
    refreshdata(p)
    pause(time)
    %t=3
    hr = 90;   bf = 0.2991;
        % bf = a*hr^2 + b*hr + c
    c = bf -  a*hr^2 - b*hr;
	bf_plot = a*hr_plot.^2+b*hr_plot+c;
    refreshdata(p)
    pause(time)
a
b
c
end
% repeat this many (inf) times and it will lead to the same result.
