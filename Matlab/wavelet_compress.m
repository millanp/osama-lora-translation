clc
close all
clear all

fi_1 = fopen('SF8_new/5_tx');
x_inter_1 = fread(fi_1, 'float32');
fclose(fi_1);
x_1_real = x_inter_1(1:2:end);
x_1_imag =  x_inter_1(2:2:end);
x_2 = x_1_real + 1i*x_1_imag;
wavelet_name = 'db5';

[thr,sorh,keepapp] = ddencmp('cmp','wv',x_1_real);
[XC_r,CXC_r,LXC_r,PERF0_r,PERFL2_r] = wdencmp('gbl',x_1_real,wavelet_name,7,thr*100,sorh,keepapp);

[thr,sorh,keepapp] = ddencmp('cmp','wv',x_1_imag);
[XC_i,CXC_i,LXC_i,PERF0_i,PERFL2_i] = wdencmp('gbl',x_1_imag,wavelet_name,7,thr*100,sorh,keepapp);

y_1_real = waverec(CXC_r,LXC_r,wavelet_name);
y_1_imag = waverec(CXC_i,LXC_i,wavelet_name);
y = y_1_real+ 1i*y_1_imag;