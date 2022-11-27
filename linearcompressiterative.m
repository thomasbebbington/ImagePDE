im = imread("image.bmp");
gray = rgb2gray(im);
gray = imgaussfilt(gray);

global h w savepixelscount indexes lambda tau;

w = size(im,2);
h = size(im,1);
tau = 0.2;

compressionratio = 0.2;

totalpixels = w * h;
savepixelscount = ceil(totalpixels * compressionratio);

randomindices = randsample(((h-2)*(w-2)),savepixelscount);
[Y,X] = ind2sub([h-2;w-2], randomindices);
Y = Y + 1;
X = X + 1;
indexes = cat(2,X,Y);

u = zeros(1,h*w);

for k = 1:savepixelscount
    u(sub2ind([h;w],indexes(k,2),indexes(k,1))) = gray(indexes(k,2),indexes(k,1));
end

avg = mean(gray,'all');
avgsaved = mean(u,'all');

u0 = u;
uncompressed = zeros(h,w);
A = generatematrix;
for count = 1:1000
   
    nu = transpose((speye(h*w)+tau*A)*transpose(u));
    u = nu;
    for k = 1:savepixelscount
        u(sub2ind([h;w],indexes(k,2),indexes(k,1))) = gray(indexes(k,2),indexes(k,1));
    end
    for i = 1:(w)
        for j = 1:(h)
            uncompressed(j,i) = u(sub2ind([h;w],j,i));
        end
    end
    uncompressed = cast(uncompressed,'uint8');

    imwrite(uncompressed,strcat("images/",string(count),".bmp"));

    disp(count);
end

% for k = 1:savepixelscount
%     u(sub2ind([h;w],indexes(k,2),indexes(k,1))) = gray(indexes(k,2),indexes(k,1));
% end

for i = 1:(w)
    for j = 1:(h)
        uncompressed(j,i) = u(sub2ind([h;w],j,i));
    end
end

selected = zeros(h,w);
for i = 1:(w)
    for j = 1:(h)
        selected(j,i) = u0(sub2ind([h;w],j,i));
    end
end

selected = cast(selected,'uint8');
imshow(gray);
figure, imshow(uncompressed);
% figure, imshow(selected);

function A = generatematrix
    global h w savepixelscount indexes tau;

    is = zeros(1,h*w);
    js = zeros(1,h*w);
    vs = zeros(1,h*w);
    index = 1;

    for i = 2:h-1
        for j = 2:w-1
            for deltax = -1:2:1
                is(index) = sub2ind([h;w],i,j+deltax);
                js(index) = sub2ind([h;w],i,j);
                vs(index) = 1;
                index = index + 1;
            end
            for deltay = -1:2:1
                is(index) = sub2ind([h;w],i+deltay,j);
                js(index) = sub2ind([h;w],i,j);
                vs(index) = 1;
                index = index + 1;
            end
            is(index) = sub2ind([h;w],i,j);
            js(index) = sub2ind([h;w],i,j);
            vs(index) = -4;
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

