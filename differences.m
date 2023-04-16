im = imread("peppers_gray.tif");
%gray = rgb2gray(im);
gray = im(:,:,1);

Gx = [-1 0 1; -2 0 2; -1 0 1];
Gy = [1 2 1; 0 0 0; -1 -2 -1];

imx = conv2(gray,Gx,'valid');
imy = conv2(gray,Gy,'valid');

J = (imx.^2 + imy.^2).^(1/2);
J = (255/max(J,[],'all'))*J;
J = cast(J,'uint8');

imshow(gray);
figure, imshow(J);

grayblur = imgaussfilt(gray,1);

imblurx = conv2(grayblur,Gx,'valid');
imblury = conv2(grayblur,Gy,'valid');

Jblur = (imblurx.^2 + imblury.^2).^(1/2);
Jblur = (255/max(Jblur,[],'all'))*Jblur;
Jblur = cast(Jblur,'uint8');

figure, imshow(Jblur);