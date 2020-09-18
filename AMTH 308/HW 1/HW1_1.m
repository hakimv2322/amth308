% This program does a 1-D Haar Transform.
close all
clear all

% . means apply to each element
% Build the s vector
s(1:7) = 70 - [1:7].^2;
s(8:16) = 70 - [1:9].^2;

% Encode
v = haar1D_encode(s,1);

% Plots
figure(1);plot([1:16],s,'k-+');
ylim([-40 100])
title('Input Signal')

figure(2);plot([1:16],v,'ro');
ylim([-40 100])
title('1D Haar Transform, Level 1')


function v = haar1D_encode(s, level)
% The function assumes length(s) is a power of 2.
% Also, s must be a row vector.
len = length(s);
% Build filter matrix
Q=[1 1;1 -1];
I = eye(len);
H = kron(I(1:len/2,1:len/2),Q)/sqrt(2);
% Build permutation matrix
PT = I([1:2:len],:);
PB = I([2:2:len],:);
% Encode input vector
u = s'; % ' means transpose
dummy = min(level, log2(len));
for j = 1:dummy
    P  = [PT(1:len/2,1:len);PB(1:len/2,1:len)];
    H  = H(1:len,1:len);
    u(1:len)=P*H*u(1:len);
    len = len/2;
end
v=u';
end
