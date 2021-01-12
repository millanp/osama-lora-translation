clear all;
clc;
close all;
rng(1)
N = 5;
k = 10;
lamda_n = 10;
x = rand(N,lamda_n*50);
T = -log(x) / lamda_n;
mean(T,2)
tot_time = sum(T,2) + 40e-3*size(T,2)