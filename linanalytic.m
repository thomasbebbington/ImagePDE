im = imread("image.bmp");
gray = rgb2gray(im);

global h w savepixelscount indexes;

w = size(im,2);
h = size(im,1);

compressionratio = 0.2;

totalpixels = w * h;
savepixelscount = ceil(totalpixels * compressionratio);

randomindices = randsample(((h-2)*(w-2)),savepixelscount);
[Y,X] = ind2sub([h-2;w-2], randomindices);
Y = Y + 1;
X = X + 1;
indexes = cat(2,X,Y);

avg = mean(gray,'all');

u = zeros(h,w);
for k = 1:savepixelscount
    u(indexes(k,2),indexes(k,1)) = gray(indexes(k,2),indexes(k,1));
end

for i = 1:100000
    u = imgaussfilt(u,3,Padding="symmetric");
    disp(i);
end

for k = 1:savepixelscount
    u(indexes(k,2),indexes(k,1)) = gray(indexes(k,2),indexes(k,1));
end

uncompressed = zeros(h,w);

for i = 1:(w)
    for j = 1:(h)
        uncompressed(j,i) = u(j,i);
    end
end

uncompressed = cast(uncompressed,'uint8');
imwrite(uncompressed,"linear.bmp");
imshow(gray);
figure, imshow(uncompressed);
