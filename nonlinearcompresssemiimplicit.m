im = imread("image.bmp");
gray = rgb2gray(im);
%gray = im;
global h sd w savepixelscount indexes lambda tau savedpixels;

sd = 0.5;
w = size(im,2);
h = size(im,1);
lambda = 0.1;
tau = 10000;

compressionratio = 0.2;

totalpixels = w * h;
savepixelscount = ceil(totalpixels * compressionratio);

randomindices = randsample(((h-2)*(w-2)),savepixelscount);
[Y,X] = ind2sub([h-2;w-2], randomindices);
Y = Y + 1;
X = X + 1;

savedpixels = false(h,w);


indexes = cat(2,X,Y);

u = zeros(1,h*w);
%u = u + mean(gray,'all');

for k = 1:savepixelscount
    u(sub2ind([h;w],indexes(k,2),indexes(k,1))) = gray(indexes(k,2),indexes(k,1));
    savedpixels(indexes(k,2),indexes(k,1)) = 1;
end

u = u./256;

u0 = u;
uncompressed = zeros(h,w);

count = 1;
changeperpixel = 10;

while count < 200 %&& changeperpixel > 0.001
    uncompressed = vectortoimage(u);
    uncompressed = uncompressed.*256;
    uncompressed = cast(uncompressed,'uint8');

    imwrite(uncompressed,strcat("images/",string(count),".bmp"));
    A = generatematrix(u);
    nu = transpose((speye(h*w)-tau*A)\transpose(u));
    
    uncompressed = vectortoimage(nu);
    uncompressed = uncompressed.*256;
    uncompressed = cast(uncompressed,'uint8');

    imwrite(uncompressed,strcat("images/",string(count),".bmp"));
    
    for k = 1:savepixelscount
        nu(sub2ind([h;w],indexes(k,2),indexes(k,1))) = u0(sub2ind([h;w],indexes(k,2),indexes(k,1)));
    end
    changeperpixel = sum(abs(nu-u),'all')/totalpixels;
    disp(changeperpixel);

    u = nu;
    
    disp(count);
    count = count + 1;
end

selected = vectortoimage(u0);

selected = cast(selected,'uint8');
imshow(gray);
figure, imshow(uncompressed);
% figure, imshow(selected);

function A = generatematrix(u)
    global h w sd savedpixels;

    is = zeros(1,h*w);
    js = zeros(1,h*w);
    vs = zeros(1,h*w);
    index = 1;

    dxu = zeros(h,w);
    dyu = zeros(h,w);

    UImage = vectortoimage(u);
    UImage = imgaussfilt(UImage,sd,Padding="symmetric");

    for i = 2:h-1
        for j = 2:w-1
            dxu(i,j) = (UImage(i,j+1) - UImage(i,j-1))/2;
            dyu(i,j) = (UImage(i+1,j) - UImage(i-1,j))/2;
        end
    end

    gu = zeros(h,w);

    for i = 2:h-1
        for j = 2:w-1
            gu(i,j) = g((dxu(i,j)^2) + (dyu(i,j)^2));
        end
    end

    for i = 2:h-1
        for j = 2:w-1
            centre = 0;
            sum = 0;
            for deltax = -1:2:1
                is(index) = sub2ind([h;w],i,j+deltax);
                js(index) = sub2ind([h;w],i,j);
                sum = (gu(i,j) + gu(i,j+deltax))/2;
                centre = centre - sum;
                vs(index) = sum;
                index = index + 1;
            end
            for deltay = -1:2:1
                is(index) = sub2ind([h;w],i+deltay,j);
                js(index) = sub2ind([h;w],i,j);
                sum = (gu(i,j) + gu(i+deltay,j))/2;
                centre = centre - sum;
                vs(index) = sum;
                index = index + 1;
            end
            is(index) = sub2ind([h;w],i,j);
            js(index) = sub2ind([h;w],i,j);
            vs(index) = centre;
            index = index + 1;
        end
    end
    
    for i = 2:h-1
        is(index) = sub2ind([h;w],i,1);
        js(index) = sub2ind([h;w],i,1);
        vs(index) = -4;
        index = index + 1;
        for deltay = -1:2:1
            is(index) = sub2ind([h;w],i+deltay,1);
            js(index) = sub2ind([h;w],i,1);
            vs(index) = 1;
            index = index + 1;
        end
        is(index) = sub2ind([h;w],i,1+1);
        js(index) = sub2ind([h;w],i,1);
        vs(index) = 2;
        index = index + 1;

        is(index) = sub2ind([h;w],i,w);
        js(index) = sub2ind([h;w],i,w);
        vs(index) = -4;
        index = index + 1;
        for deltay = -1:2:1
            is(index) = sub2ind([h;w],i+deltay,w);
            js(index) = sub2ind([h;w],i,w);
            vs(index) = 1;
            index = index + 1;
        end
        is(index) = sub2ind([h;w],i,w-1);
        js(index) = sub2ind([h;w],i,w);
        vs(index) = 2;
        index = index + 1;
    end

    for j = 2:w-1
        is(index) = sub2ind([h;w],1,j);
        js(index) = sub2ind([h;w],1,j);
        vs(index) = -4;
        index = index + 1;
        for deltax = -1:2:1
            is(index) = sub2ind([h;w],1,j+deltax);
            js(index) = sub2ind([h;w],1,j);
            vs(index) = 1;
            index = index + 1;
        end
        is(index) = sub2ind([h;w],1+1,j);
        js(index) = sub2ind([h;w],1,j);
        vs(index) = 2;
        index = index + 1;

        is(index) = sub2ind([h;w],h,j);
        js(index) = sub2ind([h;w],h,j);
        vs(index) = -4;
        index = index + 1;
        for deltax = -1:2:1
            is(index) = sub2ind([h;w],h,j+deltax);
            js(index) = sub2ind([h;w],h,j);
            vs(index) = 1;
            index = index + 1;
        end
        is(index) = sub2ind([h;w],h-1,j);
        js(index) = sub2ind([h;w],h,j);
        vs(index) = 2;
        index = index + 1;
    end

    is(index) = sub2ind([h;w],1,1);
    js(index) = sub2ind([h;w],1,1);
    vs(index) = -4;
    index = index + 1;
    is(index) = sub2ind([h;w],1,2);
    js(index) = sub2ind([h;w],1,1);
    vs(index) = 2;
    index = index + 1;
    is(index) = sub2ind([h;w],2,1);
    js(index) = sub2ind([h;w],1,1);
    vs(index) = 2;
    index = index + 1;

    is(index) = sub2ind([h;w],h,1);
    js(index) = sub2ind([h;w],h,1);
    vs(index) = -4;
    index = index + 1;
    is(index) = sub2ind([h;w],h-1,1);
    js(index) = sub2ind([h;w],h,1);
    vs(index) = 2;
    index = index + 1;
    is(index) = sub2ind([h;w],h,2);
    js(index) = sub2ind([h;w],h,1);
    vs(index) = 2;
    index = index + 1;

    is(index) = sub2ind([h;w],1,w);
    js(index) = sub2ind([h;w],1,w);
    vs(index) = -4;
    index = index + 1;
    is(index) = sub2ind([h;w],1,w-1);
    js(index) = sub2ind([h;w],1,w);
    vs(index) = 2;
    index = index + 1;
    is(index) = sub2ind([h;w],2,w);
    js(index) = sub2ind([h;w],1,w);
    vs(index) = 2;
    index = index + 1;

    is(index) = sub2ind([h;w],h,w);
    js(index) = sub2ind([h;w],h,w);
    vs(index) = -4;
    index = index + 1;
    is(index) = sub2ind([h;w],h,w-1);
    js(index) = sub2ind([h;w],w,w);
    vs(index) = 2;
    index = index + 1;
    is(index) = sub2ind([h;w],h-1,w);
    js(index) = sub2ind([h;w],h,w);
    vs(index) = 2;
    index = index + 1;
    
    A = sparse(js,is,vs,h*w,h*w);
end

function val = g(x)
    global lambda;
    val = 1/(1 + (x/(lambda^2)));
end

function img = vectortoimage(u)
    global h w;
    img = zeros(h,w);
    for i = 1:(w)
        for j = 1:(h)
            img(j,i) = u(sub2ind([h;w],j,i));
        end
    end
end

function vec = imagetovector(I)
    global h w;
    vec = zeros(h*w);
    for i = 1:(w)
        for j = 1:(h)
            vec(sub2ind([h;w],j,i)) = I(j,i);
        end
    end
end
