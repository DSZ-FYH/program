clc;
clear;
close all;
% 
% GPS = importdata('.\实验5-GPS-INS组合跑车实验\GPS数据.dat');
% GPS_data_length = length(GPS);
% GPS数据差值
% GPS = interp1(1:GPS_data_length, GPS, 1:0.005:GPS_data_length, 'spline');
% GPS(1,:) = [];
% 
% save("GPS插值数据.dat","GPS",'-ascii','-double')

A = [1+i;2-2i;3+3i;4;5;6];
A
A.'
