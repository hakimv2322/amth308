% Applies a full 2D Haar Transform, then
% sets to zero all pixels below threshold, and
% reconstructs the image.

close all
clear all

len = 256;

% Reads the original image.
A = imread('Pathfinder.jpg');
A = rgb2gray(A);
A = imresize(A,[len,len],'bicubic');
A = double(A);
Or = A;
% image(A);axis equal;colormap gray(256)

% Find the number of bytes of the original image
% s=dir('/Users/VictorHakim/Dropbox/Old Stuff/MATLAB/AMTH 308/Pathfinder.jpg');
% orig_bytes = s.bytes;
orig_bytes = len^2;

% Choose level of transform
lvl = 100;

% Choose cutoff fraction of pixels to be zero
cutoff = 0.9;

% Apply 2D Haar Transform
A = haar2D_encode(A,lvl);

% Set pixels below cutoff to be zero
[A, th] = setZero(A, cutoff);

% Quantize the image
[y,s,c] = logQuant(A(:), th, 8);

% Compress the image and sign vectors; get ratio
ratio = compress_lossless(y,s,orig_bytes);

% Uncompress the image and sign vectors
[y, s] = uncompress_lossless('Bins','Sign');

% Dequantize the image
A = deQuant(y, s, c, len);

% Apply and display inverse 2D Haar Transform
A = haar2D_decode(A,lvl);
image(A);axis equal;colormap gray(256)

% Compute the PSNR
PSNR = psnr(Or, A, len-1);


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

% Function to set pixels below cutoff to zero
function [C, th] = setZero(B, cut)
% cut should be a number between 0 and 1
[numRows,numCols] = size(B);
sz = length(B(:));
Z = sort(abs(B(:)));
thresh = Z(max(1,ceil(cut*sz)));
for it1 = 1:numRows
    for it2 = 1:numCols
        if abs(B(it1, it2)) <= thresh
            B(it1, it2) = 0;
        end
    end
end
C = B; th = thresh;
end

% Function for logarithmic quantization
function [y,s,c] = logQuant(x, th, bits)
% x is a vector, th is the threshold
NX=size(x); NP = 2^bits; % number of bins
a=abs(x); k=1;
for n=1:NX
    if a(n)>th; s(k)=sign(x(n));k=k+1;end
end
MX = max(a); c=zeros(NP,1);
p=zeros(NP-1,1);
c(1)=0.; c(NP)=MX;
p(1)=th; d=(MX/th)^(1/(NP-1));
for n=2:NP-1
    p(n)=th*d^(n-1); c(n)=p(n);
end
y = uint8(quantiz(a,p));
end

% Function to de-quantize
function B = deQuant(y, s, c, len)
BQ = c(y(:)+1);
k = 1;
s = 2*s - 1; % Include this only if gzip is used!
for n=1:len^2
    if BQ(n) ~= 0
        BQ(n)=BQ(n)*s(k);
        k = k+1;
    end
end
B = reshape(BQ,[len,len]);
end

% Lossless compression of the bins and sgn vectors
function ratio = compress_lossless(bins,sgn,original_bytes)
working_path = pwd;
%write index array (bins) to binary file 'Bins'
FILE1='Bins';fid=fopen(FILE1,'w');count=fwrite(fid,bins);status=fclose(fid);
%write array sgn to binary file  ?Sign?
FILE2='Sign';fid=fopen(FILE2,'w');count=fwrite(fid,sgn);status=fclose(fid);
%gzip both files
gzip(FILE1);gzip(FILE2);
%Size of 'Bins' and 'Sign' in bytes after gzip.
FILE1_BYTES=strcat(working_path,'/',FILE1,'.gz');
s=dir(FILE1_BYTES);compressed1_bytes = s.bytes;
FILE2_BYTES=strcat(working_path,'/',FILE2,'.gz');
s=dir(FILE2_BYTES);compressed2_bytes = s.bytes;
%Ratio between sizes of original and compressed files in bytes
ratio = original_bytes/(compressed1_bytes+compressed2_bytes);
end

% Uncompress the bins and sgn vectors
function [bins,sgn] = uncompress_lossless(FILE1,FILE2)
working_path = pwd;
GZIP1=strcat(FILE1,'.gz');GZIP2=strcat(FILE2,'.gz');
gunzip(GZIP1);gunzip(GZIP2);
fid=fopen(FILE1,'r','l');bins=fread(fid);status=fclose(fid);
fid=fopen(FILE2,'r','l');sgn=fread(fid);status=fclose(fid);
end

% Compute the PSNR
function psn = psnr(Or, R, max)
% O is original matrix
% R is reconstructed matrix
mse = 0;
[numRows,numCols] = size(Or);
for it1 = 1:numRows
    for it2 = 1:numCols
        mse = mse + (Or(it1,it2) - R(it1,it2))^2;
    end
end
mse = mse/(numRows*numCols);
psn = 10*log10((max^2)/mse);
end


