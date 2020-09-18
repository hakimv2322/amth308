% This program converts an image into 256x256 grayscale.

close all
clear all

A = imread('Pathfinder.jpg');
A = rgb2gray(A);
A = imresize(A,[256,256],'bicubic');
A = double(A);
image(A);axis equal;colormap gray(256) 
