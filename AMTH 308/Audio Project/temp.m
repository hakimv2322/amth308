clear; clc;
% load myRecording.m4a;filename='myRecording.wav';
% audiowrite(filename,y,Fs);
[y,Fs] = audioread('myRecording.m4a');
audiowrite('myRecording.wav',y,Fs);
x=y(1:32768,1);
% x=y(1:65536,1);
% x=y(1:119744,1);
sound(x,Fs)

% clear;clc;
% load handel.mat;filename='handel.wav';
% audiowrite(filename,y,8192);
% x=y(1:32768);
% sound(x,Fs)