clear all;clc;close all;
im1 = imread('2_sheep.jpg');
im2 = imread('2_grass.jpg');
im2 = imresize(im2,[864 1152]);
% figure,imshow(im2);

msk = roipoly(im1);
msk = double(msk); % Mask wrt source

% Texture
% im2 = zeros(100,100,3);
% im2(:,:,1) = 0;im2(:,:,2) = 0;im2(:,:,3) = 0;
% im2 = uint8(im2);
% msk = ones(size(im1,1),size(im1,2));
% msk(1,:) = 0;msk(:,1)=0;msk(size(msk,1),:)=0;msk(:,size(msk,2))=0;
% i1 = 51;j1 = 51;

[r,c] = find(msk > 0);
rmin = min(r);rmax = max(r);
cmin = min(c);cmax = max(c);
msk2 = msk(rmin:rmax,cmin:cmax); % Cropped mask

h=figure;
h,imshow(im2);
[j1,i1] = getpts(h);
j1 = round(j1);i1 = round(i1);

if mod(size(msk2,1),2) == 0
    r1 = size(msk2,1)/2;
    r2 = r1-1;
else
    r1 = floor(size(msk2,1)/2);
    r2 = r1;
end
if mod(size(msk2,2),2) == 0
    c1 = size(msk2,2)/2;
    c2 = c1-1;
else
    c1 = floor(size(msk2,2)/2);
    c2 = c1;
end

im1 = double(im1);
im3 = zeros(size(im2));
im4 = zeros(size(im2,1),size(im2,2));
im4(i1-r1:i1+r2,j1-c1:j1+c2) = msk2; % Mask wrt target
for c = 1:3
    temp = im1(:,:,c).*msk;
    im3(i1-r1:i1+r2,j1-c1:j1+c2,c) = temp(rmin:rmax,cmin:cmax);
end
im3 = uint8(im3); % Selected part of source pasted into selected area of target
% figure,imshow(im3);

T0 = blend(im2,im3,im4);
figure,imshow(T0);