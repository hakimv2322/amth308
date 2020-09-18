% Sample image given a one-level Haar encoding

close all
clear all

% Create the matrix
sqDim = 256; % must be even
X = eye(sqDim-2) + flip(eye(sqDim-2));
Y = [ones(sqDim-2,1) X ones(sqDim-2,1)];
A = [ones(1,sqDim); Y; ones(1,sqDim)];
A = (sqDim-1)*(-A + ones(sqDim));
A = double(A);
% image(A);axis equal;colormap gray(256)

% Apply and display one-level 2D Haar Transform
B = haar2D_encode(A, 1);
image(B);axis equal;colormap gray(256)

% 2D Haar Transform function
function B = haar2D_encode(A, level)
% The function assumes A is a square and power of 2.
len = length(A);
% Build filter matrix
Q=[1 1;1 -1];
I = eye(len);
H = kron(I(1:len/2,1:len/2),Q)/sqrt(2);
% Build permutation matrix
PT = I([1:2:len],:);
PB = I([2:2:len],:);
% Encode input vector
dummy = min(level, log2(len));
for j = 1:dummy
    P  = [PT(1:len/2,1:len);PB(1:len/2,1:len)];
    H  = H(1:len,1:len);
    A(1:len,1:len)=P*H*A(1:len,1:len)*H'*P';
    len = len/2;
end
B = A;
end

