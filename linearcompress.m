im = imread("image.bmp");
gray = rgb2gray(im);

w = size(im,2);
h = size(im,1);

n = min([w h]);
gray = gray(1:n,1:n);

compressionratio = 0.5;

savepixelscount = ceil(n^2 * compressionratio);

S = randsample(n^2,savepixelscount);
U = setdiff(1:n^2,S);

M=gallery('poisson',n);
A = M(1:n^2,1:n^2);

b = zeros(n^2,1);
b(S) = gray(S);

b = -A*b;
b = b(U);
AA = A(U,U);

tic
img = AA\b;
%img = cgs(AA,b,1e-6,200);
%img = cgs(AA,b,1e-6,200,D);
toc

uncompressed = zeros(n);

uncompressed(U) = img;
uncompressed(S) = gray(S);

uncompressed = cast(uncompressed,'uint8');
imwrite(uncompressed,"linear.bmp");
figure, imshow(gray);
figure, imshow(uncompressed);