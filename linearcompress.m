im = imread("image.bmp");
gray = rgb2gray(im);

w = size(im,2);
h = size(im,1);

n = min([w h]);
gray = gray(1:n,1:n);

compressionratio = 0.02;

savepixelscount = ceil(n^2 * compressionratio);

S = randsample(n^2,savepixelscount);
U = setdiff(1:n^2,S);

M=gallery('poisson',n);
A = M(1:n^2,1:n^2);

b = zeros(1,n^2);
b(S) = gray(S);

b = transpose(-A*transpose(b));
AA = sparse(n^2,n^2);
I = speye(n^2,n^2);
AA(:,S) = I(:,S); 
AA(:,U) = A(:,U);

ds = spdiags(AA,0);
D = spdiags(ds,0,n^2,n^2);


tic
img = AA\transpose(b);
%img = cgs(AA,transpose(b),1e-6,200);
%img = cgs(AA,transpose(b),1e-6,200,D);
toc

uncompressed = zeros(n);

uncompressed(U) = img(U);
uncompressed(S) = gray(S);

uncompressed = cast(uncompressed,'uint8');
imwrite(uncompressed,"linear.bmp");
figure, imshow(gray);
figure, imshow(uncompressed);