im = imread("image.png");
gray = rgb2gray(im);
%gray = im(:,:,1);

w = size(im,2);
h = size(im,1);

n = min([w h]);
gray = gray(1:n,1:n);

compressionratio = 0.2;

savepixelscount = ceil(n^2 * compressionratio);
unsavepixelcount = n^2 - savepixelscount;

nc = floor((n-1)/2);

S = randsample(n^2,savepixelscount);

U = setdiff(1:n^2,S);

M=gallery('poisson',n);
A = M(1:n^2,1:n^2);

b = zeros(n^2,1);
b(S) = gray(S);

ob = b;

b = -A*b;
b = b(U);

AA = A(U,U);

R = generateRestrict(nc,n);
R = R(:,U);
P = transpose(R);

AAC = R*AA*P;

uf = zeros(unsavepixelcount,1);

uncompressed = zeros(n^2,1);
uncompressed(U) = uf;
uncompressed(S) = gray(S);
uncompressed = reshape(uncompressed,[n n]);
uncompressed = uncompressed(2:n-1,2:n-1);
uncompressed = cast(uncompressed,'uint8');
figure, imshow(uncompressed);

AAL = tril(AA,0);
AAU = triu(AA,1);

AACL = tril(AAC,0);
AACU = triu(AAC,1);

tic
for iterationcount = 1:6
    for relaxationcount = 1:2
        uf = relax(uf,AAL,AAU,b);
    end
    
    rf = b - (AA*uf);
    
    rc = R*rf;
    
    ec = zeros(nc^2,1); 
    for relaxationcount = 1:3
        ec = relax(ec,AACL,AACU,rc);
    end
    
    ef = P*ec;

    uf = uf + ef;
end
for relaxationcount = 1:2
    uf = relax(uf,AAL,AAU,b);
end
toc

tic
ufcg = cgs(AA,b);
toc

uncompressed = zeros(n^2,1);
uncompressed(U) = uf;
uncompressed(S) = gray(S);
uncompressed = reshape(uncompressed,[n n]);
uncompressed = uncompressed(2:n-1,2:n-1);

uncompressedcg = zeros(n^2,1);
uncompressedcg(U) = ufcg;
uncompressedcg(S) = gray(S);
uncompressedcg = reshape(uncompressedcg,[n n]);
uncompressedcg = uncompressedcg(2:n-1,2:n-1);

uncompressedcg = cast(uncompressedcg,'uint8');
figure, imshow(uncompressedcg);

uncompressed = cast(uncompressed,'uint8');
imwrite(uncompressed,"linear.bmp");
figure, imshow(gray);
figure, imshow(uncompressed);

function R = generateRestrict(nc,nf)
    is = zeros(1,nc^2);
    js = zeros(1,nf^2);
    vs = zeros(1,nf*nc);
    index = 1;
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
        end
    end
    
    is = nonzeros(is);
    js = nonzeros(js);
    vs = nonzeros(vs);

    R = sparse(is,js,vs,nc^2,nf^2);
end

function u = relax(u,L,U,b)
    u = L\(b - U*u);
    u = transpose(L)\(b - transpose(U)*u);
end

function u = relaxj(u,L,U,D,b)
    u = D\(b-(L+U)*u);
end