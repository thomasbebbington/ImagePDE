im = imread("image.png");
gray = rgb2gray(im);

w = size(im,2);
h = size(im,1);

n = min([w h]);
gray = gray(1:n,1:n);

compressionratio = 0.3;

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

uf = zeros(unsavepixelcount,1);

v = zeros(n^2,1);

AAC = gallery('poisson',nc);
AAC = ((n^2))*AAC;

tic
for iterationcount = 1:4
    for relaxationcount = 1:2
        uf = relax(uf,AA,b);
    end
    
    rf = b - (AA*uf);

    v(U) = rf;

    rc = restrict(v,nc,n);

    ec = AAC\rc;
    
    ef = prolong(ec,n,nc);
    
    uf = uf + ef(U);
end

for i = 1:2
    uf = relax(uf,AA,b);
end
toc

uncompressed = zeros(n);

uncompressed(U) = uf;
uncompressed(S) = gray(S);

disp(norm(double(gray(2:n-1,2:n-1))-uncompressed(2:n-1,2:n-1)));

uncompressed = cast(uncompressed,'uint8');
imwrite(uncompressed,"linear.bmp");
figure, imshow(gray);
figure, imshow(uncompressed);

function uc = restrict(uf,nc,nf)
    UF = reshape(uf,[nf nf]);
    UC = zeros(nc,nc);
    for i = 1:nc
        for j = 1:nc
            sum = (UF((2*i)-1,(2*j)-1) + UF((2*i)-1,(2*j)+1) + UF((2*i)+1,(2*j)-1) + UF((2*i)+1,(2*j)+1) ...
                + 2*(UF((2*i),(2*j)-1) + UF((2*i),(2*j)+1) + UF((2*i)-1,(2*j)) + UF((2*i)+1,(2*j))) ...
                + 4*UF((2*i),(2*j)));
            UC(i,j) = (1/16)*sum;
        end
    end
    uc = reshape(UC,[nc^2 1]);
end

function uf = prolong(uc,nf,nc)
    UC = reshape(uc,[nc nc]);
    UF = zeros(nf,nf);

    % Horizontals and Verticals
    for i = 1:nc-1
        for j = 1:nc-1
            UF((2*i)+1,2*j) = (1/2)*(UC(i,j) + UC(i+1,j));
            UF((2*i),(2*j)+1) = (1/2)*(UC(i,j) + UC(i,j+1));
        end
    end

    % Diagonals
    for i = 1:nc-1
        for j = 1:nc-1
            UF((2*i)+1,(2*j)+1) = (1/4)*(UC(i,j) + UC(i+1,j+1) + UC(i+1,j) + UC(i,j+1));
        end
    end

    % Adjacent Sides
    for i = 1:nc
        UF((2*i),1) = UC(i,1);
        UF(1,(2*i)) = UC(1,i);
        UF((2*i),2*(nc)+1) = UC(i,nc);
        UF(2*(nc)+1,(2*i)) = UC(nc,i);
    end

    % Non-Adjacent Sides
    for i = 1:nc-1
        UF((2*i)+1,1) = (1/2)*(UC(i,1) + UC(i+1,1));
        UF((2*i)+1,2*(nc)+1) = (1/2)*(UC(i,nc) + UC(i+1,nc));

        UF(1,(2*i)+1) = (1/2)*(UC(1,i) + UC(1,i+1));
        UF(2*(nc)+1,(2*i)+1) = (1/2)*(UC(nc,i) + UC(nc,i));
    end

    % Corners
    UF(1,1) = UC(1,1);
    UF((2*nc)+1,1) = UC(nc,1);
    UF(1,(2*nc)+1) = UC(1,nc);
    UF((2*nc)+1,(2*nc)+1) = UC(nc,nc);

    %Additional
    if(2*nc ~= nf)
        for i = 1:nf
            UF(i,nf) = UF(i,nf-1);
            UF(nf,i) = UF(nf-1,i);
        end
        UF(nf,nf) = UF(nf-1,nf-1);
    end

    uf = reshape(UF,[nf^2 1]);
end

function u = relax(u,A,b)
    L = tril(A,0);
    U = triu(A,1);

    u = L\(b - U*u);
    u = (transpose(L))\(b - ((transpose(U)*(u))));
end