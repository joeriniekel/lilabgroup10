%%
clc
% 
% amp=0.5;
% duration=0.19;
% fs=20500  % sampling frequency
% freq=800;
% values=0:1/fs:duration;
% a=amp*sin(2*pi* freq*values);
% sound(a);
% 
% pause(2)
% 
% freq=1200;
% values=0:1/fs:duration;
% a=amp*sin(2*pi* freq*values);
% sound(a);

%%
dt = 0.2
cycle_time = 0;
  amp = 0.5;
  duration = 0.2;
  Fs = 8192;  % sampling frequency
  values=0:1/Fs:duration;
  len = length(values)
  decreasing_amp = amp:-1*amp/len:0;
  
      f = 800 + (cycle_time*1).^2;
      a2 = sin(2*pi* f * values);
      if length(decreasing_amp) > length(values), disp('too long');end;
       a2 = a2 .* decreasing_amp(1:len);
      sound(a2,Fs);
      
      
      