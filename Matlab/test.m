% load leleccum; indx = 2600:3100;
% x = leleccum(indx);
clear all
clc
X = [1:100] + i*[100:-1:1]
x = awgn(X,0);
x_r = real(x);
x_i = imag(x);
%%
% plot(real(x))
[thr,sorh,keepapp] = ddencmp('den','wv',x);
xd = wdencmp('gbl',x,'db3',2,thr,sorh,keepapp);
subplot(211)
plot(x); title('Original Signal');
subplot(212)
plot(xd); title('Denoised Signal');
%%
clear all
close all;
load sinsin
Y = X+18*randn(size(X));
[thr,sorh,keepapp] = ddencmp('den','wv',Y);
xd = wdencmp('gbl',Y,'sym4',2,thr,sorh,keepapp);
subplot(2,2,1)
imagesc(X)
title('Original Image')
subplot(2,2,2)
imagesc(Y)
title('Noisy Image')
subplot(2,2,3)
imagesc(xd)
title('Denoised Image')