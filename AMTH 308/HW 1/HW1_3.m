% Displays a grayscale image, applies a full 2D Haar Transform,
% then applies the inverse 2D Haar Transform

close all
clear all

% Reads and displays the original image.
A = imread('Pathfinder.jpg');
A = rgb2gray(A);
A = imresize(A,[256,256],'bicubic');
A = double(A);
image(A);axis equal;colormap gray(256)

% Choose level of transform
lvl = 100;

% Apply and display 2D Haar Transform
B = haar2D_encode(A,lvl);
image(B);axis equal;colormap gray(256)

% Apply and display inverse 2D Haar Transform
C = haar2D_decode(B,lvl);
image(C);axis equal;colormap gray(256)


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

% 2D Inverse Haar Transform function
function A = haar2D_decode(B, level)
% The function assumes B is a square and power of 2.
len = length(B);
% Build filter matrix
Q=[1 1;1 -1];
I = eye(len);
H = kron(I(1:len/2,1:len/2),Q)/sqrt(2);
% Build permutation matrix
PT = I([1:2:len],:);
PB = I([2:2:len],:);
% Encode input vector
dummy = min(level, log2(len));
len = len/(2^(dummy-1));
for j = 1:dummy
    P  = [PT(1:len/2,1:len);PB(1:len/2,1:len)];
    T  = H(1:len,1:len);
    B(1:len,1:len)=T'*P'*B(1:len,1:len)*P*T;
    len = len*2;
end
A = B;
end