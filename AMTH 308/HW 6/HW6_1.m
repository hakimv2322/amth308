% This program does a 1-D Daubechies-4 Transform.
close all
clear all

% Build the time vector t
t = 1/16:1/8:32-1/16;

% Build the signal vector s
syms z
f1 = 16*z; f2 = (z-32)^2;
for it = 1:128
    z = t(it);
    s(it) = subs(f1);
end
for it = 129:256
    z = t(it);
    s(it) = subs(f2);
end

% Encode
len = length(s);
v = DaubFour_encode(s, 100);

% Choose cutoff fraction of pixels to be zero
cutoff = 0.8;

% Set pixels below cutoff to zero
v2 = setZero(v, cutoff);

% Decode the signal
s2 = DaubFour_decode(v2, 100);

% Plots
figure(1);plot(t,s,'k.');
xlim([0 32])
title('Input Signal')

figure(2);plot(t,v,'k.');
xlim([0 32])
title('1D Daubechies-4 Transform, Full Level')

figure(3);plot(t,v2,'k.');
xlim([0 32])
title('1D Daubechies-4 Transform, Full Level, 80% Cutoff')

figure(4);plot(t,s2,'k.');
xlim([0 32])
ylim([0 300])
title('Reconstructed Signal')

function T = DaubFour_Matrix(len)
% The function assumes len is a power of 2.
I = eye(len);
% Write the filter coefficients
h0 = (1/(4*sqrt(2)))*(1+sqrt(3));
h1 = (1/(4*sqrt(2)))*(3+sqrt(3));
h2 = (1/(4*sqrt(2)))*(3-sqrt(3));
h3 = (1/(4*sqrt(2)))*(1-sqrt(3));
% Build filter matrix
Q1=[h0 h1; h3 -h2]; Q2 = [h2 h3; h1 -h0];
T1 = kron(I(1:len/2,1:len/2),Q1);
offDiag = diag(ones(1,len/2-1), 1);
T2 = kron(offDiag, Q2);
corner = zeros(len/2); corner(len/2, 1) = 1;
T3 = kron(corner, Q2);
T = T1 + T2 + T3;
end

% Encode function
function v = DaubFour_encode(s, level)
% The function assumes length(s) is a power of 2.
% Also, s must be a row vector.
len = length(s);
I = eye(len);
% Build permutation matrix
PT = I([1:2:len],:);
PB = I([2:2:len],:);
% Encode input vector
u = s'; % ' means transpose
dummy = min(level, log2(len));
for j = 1:dummy
    P  = [PT(1:len/2,1:len);PB(1:len/2,1:len)];
    T  = DaubFour_Matrix(len);
    u(1:len)=P*T*u(1:len);
    len = len/2;
end
v=u';
end

% Decode function
function A = DaubFour_decode(v, level)
% The function assumes s is a row vector with length a power of 2.
len = length(v);
I = eye(len);
% Build permutation matrix
PT = I([1:2:len],:);
PB = I([2:2:len],:);
% Encode input vector
dummy = min(level, log2(len));
len = len/(2^(dummy-1));
u = v';
for j = 1:dummy
    P  = [PT(1:len/2,1:len);PB(1:len/2,1:len)];
    U  = DaubFour_Matrix(len);
    u(1:len)=U'*P'*u(1:len);
    len = len*2;
end
A = u';
end

% Function to set pixels below cutoff to zero
function C = setZero(s, cut)
% cut should be a number between 0 and 1
% s should be a vector
len = length(s);
Z = sort(abs(s));
thresh = Z(max(1,ceil(cut*len)));
for it = 1:len
        if abs(s(it)) <= thresh
            s(it) = 0;
    end
end
C = s;
end
