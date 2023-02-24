im = imread("image.bmp");
gray = rgb2gray(im);

fwx = [0 -1 1];
dxf = abs(conv2(gray,fwx,'valid'));
dxfim = cast(dxf,'uint8');
imshow(dxfim);

bwx = [-1 1 0];
dxb = abs(conv2(gray,bwx,'valid'));
dxbim = cast(dxb,'uint8');
figure, imshow(dxbim);

cx = [-1 0 1];
dxc = abs(conv2(gray,cx,'valid'));
dxcim = cast(dxc,'uint8');
figure, imshow(dxcim);

fwy = [0; -1; 1];
dyf = abs(conv2(gray,fwy,'valid'));
dyfim = cast(dyf,'uint8');
figure, imshow(dyfim);