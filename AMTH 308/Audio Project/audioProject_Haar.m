% This program takes an audio recording, converts the file into
% an array of length 32768, compresses the array using the Haar transform
% with a specified cutoff, then reconstructs the array.
% The original array's sound is played, followed immediately with
% the reconstructed array's sound.

close all
clear all

% Write the audio file into a length 32768 array
[y,Fs] = audioread('myRecording.m4a');
audiowrite('myRecording.wav',y,Fs);
s=y(1:32768,1);

% Choose cutoff fraction of pixels to be zero
cutoff = 0.94;

s = s';

% Break up the array into pieces
q = 4096;
sa = s(1:q);
sb = s(q+1:2*q);
sc = s(2*q+1:3*q);
sd = s(3*q+1:4*q);
se = s(4*q+1:5*q);
sf = s(5*q+1:6*q);
sg = s(6*q+1:7*q);
sh = s(7*q+1:8*q);

% Encode
len = length(sa);
va = Haar_encode(sa, 100);
vb = Haar_encode(sb, 100);
vc = Haar_encode(sc, 100);
vd = Haar_encode(sd, 100);
ve = Haar_encode(se, 100);
vf = Haar_encode(sf, 100);
vg = Haar_encode(sg, 100);
vh = Haar_encode(sh, 100);

% Set pixels below cutoff to zero
va2 = setZero(va, cutoff);
vb2 = setZero(vb, cutoff);
vc2 = setZero(vc, cutoff);
vd2 = setZero(vd, cutoff);
ve2 = setZero(ve, cutoff);
vf2 = setZero(vf, cutoff);
vg2 = setZero(vg, cutoff);
vh2 = setZero(vh, cutoff);


% Decode the signal
sa2 = Haar_decode(va2, 100);
sb2 = Haar_decode(vb2, 100);
sc2 = Haar_decode(vc2, 100);
sd2 = Haar_decode(vd2, 100);
se2 = Haar_decode(ve2, 100);
sf2 = Haar_decode(vf2, 100);
sg2 = Haar_decode(vg2, 100);
sh2 = Haar_decode(vh2, 100);

% Put the pieces back together
s2 = [sa2,sb2,sc2,sd2,se2,sf2,sg2,sh2];

s2 = s2';
s = s';

sound(s,Fs)
pause(1.2)
sound(s2,Fs)


function T = Haar_Matrix(len)
% The function assumes len is a power of 2.
I = eye(len);
% Build filter matrix
Q=[1 1;1 -1];
I = eye(len);
T = kron(I(1:len/2,1:len/2),Q)/sqrt(2);
end

% Encode function
function v = Haar_encode(s, level)
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
    T  = Haar_Matrix(len);
    u(1:len)=P*T*u(1:len);
    len = len/2;
end
v=u';
end

% Decode function
function A = Haar_decode(v, level)
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
    U  = Haar_Matrix(len);
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

