a = 4.2857e-04
b = -0.0464
c = 2.1500

a =  -6.6667e-05
b =   0.0230
c =  -0.4833

a = 1.0000e-05
b = 0.0040
c = 0.3000

lhr = 40

hr = -100:100;
plot(hr, a*(hr-lhr).^2 + b*(hr-lhr) + c)

plot(hr, 0.00001*(hr-lhr).^2 + 0.004*(hr-lhr) + 0.3)




amp=0.5;
fs=20500  % sampling frequency
duration=0.1;
freq=800;
values=0:1/fs:duration;
a=amp*sin(2*pi* freq*values);
sound(a);

amp=0.5;
fs=20500  % sampling frequency
duration=0.1;
freq=400;
values=0:1/fs:duration;
a=amp*sin(2*pi* freq*values);
sound(a);
