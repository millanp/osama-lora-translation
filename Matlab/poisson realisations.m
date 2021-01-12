clear all;
clc;
close all;
x = rand(1,10000);
N = 5;
k = 10;
lamda_n = 2;%k/N;
T = - log(x) / lamda_n
mean(T)