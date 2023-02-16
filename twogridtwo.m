addpath bttc
im = imread("spongebob.png");
gray = rgb2gray(im);

w = size(im,2);
h = size(im,1);

n = min([w h]);
n = 513;
gray = gray(1:n,1:n);

compressionratio = 0.1;

savepixelscount = ceil(n^2 * compressionratio);
unsavepixelcount = n^2 - savepixelscount;

nc = floor((n-1)/2);

S = randsample(n^2,savepixelscount);

tri = 60;
[Vi,Fi]=bttc_m(double(gray),tri);
ind_bttc = n.*(Vi(:,1)-1)+Vi(:,2);

edge = zeros(n,n);
edge(ind_bttc) = Vi(:,3);
figure, imshow(edge);

S = union(ind_bttc,S);

U = setdiff(1:n^2,S);

M=gallery('poisson',n);
A = M(1:n^2,1:n^2);

b = zeros(n^2,1);
b(S) = gray(S);

ob = b;

b = -A*b;
b(S) = gray(S);


AA = speye(n^2,n^2);
AA(U,U) = A(U,U);

R = generateRestrict(nc,n);
P = transpose(R);

AAC = R*AA*P;

uf = zeros(n^2,1);
uf(S) = gray(S);

uncompressed = uf;
uncompressed = reshape(uncompressed,[n n]);
uncompressed = uncompressed(2:n-1,2:n-1);
uncompressed = cast(uncompressed,'uint8');
figure, imshow(uncompressed);

tic
for iterationcount = 1:2
    for relaxationcount = 1:1
        uf = relax(uf,AA,b);
    end
    
    uS = uf(S);
    
    rf = b - (AA*uf);

    rS = rf(S);
    
    rc = R*rf;

    ec = AAC\rc;
    
    ef = P*ec;

    ufim = uf;
    ufim = reshape(ufim,[n n]);
    ufim = cast(ufim,'uint8');
    error = ef;
    error = reshape(error,[n n]);
    error = cast(error,'uint8');
    %figure, imshow([ufim error]);

    uf = uf + ef;
end
for relaxationcount = 1:5
    uf = relax(uf,AA,b);
end
toc
disp(size(S)/(n^2));

uncompressed = uf;
uncompressed = reshape(uncompressed,[n n]);
uncompressed = uncompressed(2:n-1,2:n-1);

uncompressed = cast(uncompressed,'uint8');
imwrite(uncompressed,"linear.bmp");
figure, imshow(gray);
figure, imshow(uncompressed);

function R = generateRestrict(nc,nf)
    is = zeros(1,nc^2);
    js = zeros(1,nf^2);
    vs = zeros(1,nf*nc);
    index = 1;

    %R = sparse(nc^2,nf^2);
%     row = zeros(1, 2*nf + 3);
% 
%     row(1:3) = [1 2 1];
%     row(nf+1:nf+3) = [2 4 2];
%     row((2*nf + 1):(2*nf + 3)) = [1 2 1];

    for i = 1:nc
        for j = 1:nc
            is(index) = (i-1)*nc + j;
            js(index) = (2*(i-1)*nf + 2*(j-1) + 1);
            vs(index) = 1;
            index = index + 1;

            is(index) = (i-1)*nc + j;
            js(index) = (2*(i-1)*nf + 2*(j-1) + 2);
            vs(index) = 2;
            index = index + 1;

            is(index) = (i-1)*nc + j;
            js(index) = (2*(i-1)*nf + 2*(j-1) + 3);
            vs(index) = 1;
            index = index + 1;

            is(index) = (i-1)*nc + j;
            js(index) = (2*(i-1)*nf + 2*(j-1) + nf+1);
            vs(index) = 2;
            index = index + 1;
            
            is(index) = (i-1)*nc + j;
            js(index) = (2*(i-1)*nf + 2*(j-1) + nf+2);
            vs(index) = 4;
            index = index + 1;
            
            is(index) = (i-1)*nc + j;
            js(index) = (2*(i-1)*nf + 2*(j-1) + nf+3);
            vs(index) = 2;
            index = index + 1;

            is(index) = (i-1)*nc + j;
            js(index) = (2*(i-1)*nf + 2*(j-1) + 2*nf + 1);
            vs(index) = 1;
            index = index + 1;
            
            is(index) = (i-1)*nc + j;
            js(index) = (2*(i-1)*nf + 2*(j-1) + 2*nf + 2);
            vs(index) = 2;
            index = index + 1;
            
            is(index) = (i-1)*nc + j;
            js(index) = (2*(i-1)*nf + 2*(j-1) + 2*nf + 3);
            vs(index) = 1;
            index = index + 1;

            % R(((i-1)*nc + j), ((2*(i-1)*nf + 2*(j-1) + 1):(2*(i-1)*nf + 2*(j-1) + 2*nf + 3))) = row;
        end
    end
    
    is = nonzeros(is);
    js = nonzeros(js);
    vs = nonzeros(vs);

    R = sparse(is,js,vs,nc^2,nf^2);
    R = (1/16)*R;
end


function u = relax(u,A,b)
    L = tril(A,0);
    U = triu(A,1);

    u = L\(b - U*u);
    u = (transpose(L))\(b - ((transpose(U)*(u))));
end