im = imread("imagesmall.bmp");
gray = rgb2gray(im);

w = size(im,2);
h = size(im,1);

compressionratio = 0.1;

totalpixels = w * h;
savepixelscount = ceil(totalpixels * compressionratio);

randomindices = randsample(((h-2)*(w-2)),savepixelscount);
[Y,X] = ind2sub([h-2;w-2], randomindices);
Y = Y + 1;
X = X + 1;
indexes = cat(2,X,Y);

b = zeros(1,h*w);

is = zeros(1,h*w);
js = zeros(1,h*w);
vs = zeros(1,h*w);
index = 1;

for k = 1:savepixelscount
    b(sub2ind([h;w],indexes(k,2),indexes(k,1))) = gray(indexes(k,2),indexes(k,1));
    %A(sub2ind([h;w],Y(k),X(k)),sub2ind([h;w],Y(k),X(k))) = 1;
    is(index) = sub2ind([h;w],Y(k),X(k));
    js(index) = sub2ind([h;w],Y(k),X(k));
    vs(index) = 1;
    index = index + 1;
end

lambda = 1;
g = zeros(h,w);

for i = 2:h-1
    for j = 2:w-1
        dxu = (gray(i+1,j)-gray(i-1,j))/2;
        dyu = (gray(i,j+1)-gray(i,j-1))/2;
        g(i,j) = (1/(1+((sqrt(double(dxu^2 + dyu^2)))/lambda^2)));
    end
end



for i = 1:w
    for j = 1:h
        if (j==1) || (i==1) || (j==h) || (i==w)
            is(index) = sub2ind([h;w],j,i);
            js(index) = sub2ind([h;w],j,i);
            vs(index) = 1;
            index = index + 1;
            %A(sub2ind([h;w],j,i),sub2ind([h;w],j,i)) = 1;
            b(sub2ind([h;w],j,i)) = 255;
        elseif ~ismember([i j],indexes,'rows')
            c01 = 0;
            c10 = 0;
            c11 = 0;
            c12 = 0;
            c21 = 0;

            c21 = c21 + ((g(j,i+1) + g(j,i))/2);
            c11 = c11 - ((g(j,i+1) + g(j,i))/2);
    
            c11 = c11 - ((g(j,i) + g(j,i-1))/2);
            c01 = c01 + ((g(j,i) + g(j,i-1))/2);
    
            c12 = c12 + ((g(j+1,i) + g(j,i))/2);
            c11 = c11 - ((g(j+1,i) + g(j,i))/2);
    
            c11 = c11 - ((g(j,i) + g(j-1,i))/2);
            c10 = c10 + ((g(j,i) + g(j-1,i))/2);

            %A(sub2ind([h;w],j,i),sub2ind([h;w],j-1,i)) = c10;
            js(index) = sub2ind([h;w],j,i);
            is(index) = sub2ind([h;w],j-1,i);
            vs(index) = c10;
            index = index + 1;

            %A(sub2ind([h;w],j,i),sub2ind([h;w],j+1,i)) = c12;
            js(index) = sub2ind([h;w],j,i);
            is(index) = sub2ind([h;w],j+1,i);
            vs(index) = c12;
            index = index + 1;

            %A(sub2ind([h;w],j,i),sub2ind([h;w],j,i)) = c11;
            js(index) = sub2ind([h;w],j,i);
            is(index) = sub2ind([h;w],j,i);
            vs(index) = c11;
            index = index + 1;
            if c11 == 0
                disp(i);
            end

            %A(sub2ind([h;w],j,i),sub2ind([h;w],j,i-1)) = c01;
            js(index) = sub2ind([h;w],j,i);
            is(index) = sub2ind([h;w],j,i-1);
            vs(index) = c01;
            index = index + 1;

            %A(sub2ind([h;w],j,i),sub2ind([h;w],j,i+1)) = c21;
            js(index) = sub2ind([h;w],j,i);
            is(index) = sub2ind([h;w],j,i+1);
            vs(index) = c21;
            index = index + 1;
        end
    end
end


is = is(1:index-1);
js = js(1:index-1);
vs = vs(1:index-1);
A = sparse(js,is,vs);

% for i = 1:(w)
%     for j = 1:(h)
%         if (j==1) || (i==1) || (j==h) || (i==w)
%             A(sub2ind([h;w],j,i),sub2ind([h;w],j,i)) = 1;
%             b(sub2ind([h;w],j,i)) = 255;
%         elseif A(sub2ind([h;w],j,i),sub2ind([h;w],j,i)) ~= 1
%             A(sub2ind([h;w],j,i),sub2ind([h;w],j-1,i)) = 1;
%             A(sub2ind([h;w],j,i),sub2ind([h;w],j+1,i)) = 1;
%             A(sub2ind([h;w],j,i),sub2ind([h;w],j,i)) = -4;
%             A(sub2ind([h;w],j,i),sub2ind([h;w],j,i-1)) = 1;
%             A(sub2ind([h;w],j,i),sub2ind([h;w],j,i+1)) = 1;
%         end 
%     end
% end



% Gaussian Elimination

% for r = 1:(h*w)-1
%     for r2 = r+1:(h*w)
%         ratio = A(r2,r) / A(r,r);
%         for e = r:h*w
%             A(r2,e) = A(r2,e) - ratio*A(r,e);
%         end
%         b(r2) = b(r2) - ratio*b(r);
%     end
%     disp(r/(h*w));
% end
% 
% for currentrow = (h*w)-1:-1:2
%     for nextrow = currentrow-1:-1:1
%         ratio = A(nextrow,currentrow) / A(currentrow,currentrow);
%         for column = h*w:-1:nextrow+1
%             A(nextrow,column) = A(nextrow,column) - ratio*(A(currentrow,column));
%         end
%         b(nextrow) = b(nextrow) - ratio*b(currentrow);
%     end
% end

img = A\transpose(b);
uncompressed = zeros(h,w);



% for i = 1:(w)
%     for j = 1:(h)
%         uncompressed(i,j) = b(sub2ind([h;w],j,i)) / A(sub2ind([h;w],j,i),sub2ind([h;w],j,i));
%     end
% end

for i = 1:(w)
    for j = 1:(h)
        uncompressed(j,i) = img(sub2ind([h;w],j,i));
    end
end

uncompressed = cast(uncompressed,'uint8');
imshow(gray);
figure, imshow(uncompressed);