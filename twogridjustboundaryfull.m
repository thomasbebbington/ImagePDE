nf = 101;
nc = ((nf-1)/2);

AF = full(gallery("poisson",nf));

R = generateRestrict(nc,nf);
P = transpose(R);

bFull = ones(nf+2);
AFull = gallery("poisson",nf+2);

bFull(1:nf+2,1) = 0;
bFull(:,nf+2) = 255;
bFull(1,:) = 0;
bFull(nf+2,:) = 255;

bFull = reshape(bFull, [(nf+2)^2 1]);

b = -(AFull * bFull);

b = reshape(b, [nf+2 nf+2]);
b = b(2:nf+1,2:nf+1);

b = reshape(b, [nf^2 1]);

uf = ones(nf^2,1)*128;

AC = R*AF*P;

tic
for iterationcount = 1:2
    for relaxationcount = 1:2
        uf = relax(uf,AF,b);
    end
    
    rf = b - (AF*uf);

    rc = R*rf;

    ec = AC\rc;
    
    ef = P*ec;
    
    uf = uf + ef;
end

for i = 1:2
    uf = relax(uf,AF,b);
end
toc

u = reshape(uf,[nf nf]);

u = cast(u,"uint8");

imshow(u);

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