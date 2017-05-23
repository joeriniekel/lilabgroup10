%%
clc

amp=1;
duration=0.19;
fs=20500  % sampling frequency
freq=800;
values=0:1/fs:duration;
a=amp*sin(2*pi* freq*values);
sound(a);

pause(2)

freq=1200;
values=0:1/fs:duration;
a=amp*sin(2*pi* freq*values);
sound(a);


