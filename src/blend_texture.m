clear all;clc;close all;
im1 = imread('text.jpg');
im2 = imread('desk2.jpg');

im1 = imresize(im1,0.3);
% im2 = imresize(im2,0.3);
% im2 = imresize(im2,[size(im1,1),size(im1,2)]);
% im1 = imresize(im1,[size(im2,1),size(im2,2)]);

msk = roipoly(im1);
msk = double(msk);
grey = rgb2gray(im1);
edg = edge(grey,'canny');
edg = double(edg);

[r,c] = find(msk > 0);
rmin = min(r);rmax = max(r);
cmin = min(c);cmax = max(c);
msk2 = msk(rmin:rmax,cmin:cmax); % Getting the bounding rectangle for mask

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
im3_e = zeros(size(im2,1),size(im2,2));
im4(i1-r1:i1+r2,j1-c1:j1+c2) = msk2;
temp2 = edg(:,:).*msk;
im3_e(i1-r1:i1+r2,j1-c1:j1+c2) = temp2(rmin:rmax,cmin:cmax);
im3_e = logical(im3_e);
for c = 1:3
    temp1 = im1(:,:,c).*msk;
    im3(i1-r1:i1+r2,j1-c1:j1+c2,c) = temp1(rmin:rmax,cmin:cmax); %Getting source part wherever mask is 1
end
im3 = uint8(im3);
%figure,imshow(im3);

%p2 = [x1';y1'];

%% Blend
S0 = im3;
T0 = im2;

mask0 = logical(im4);

[x,y] = find(mask0);
minx = min(x) - 1;
miny = min(y) - 1;
maxx = max(x) + 1;
maxy = max(y) + 1;
n = length(x);

mask = mask0(minx:maxx,miny:maxy);
mask_edge = im3_e(minx:maxx,miny:maxy);
mask_b = mask - imerode(mask,[0 1 0;1 1 1;0 1 0]);
mask_b=logical(mask_b);
loc = find(mask(:));
grid = zeros(size(mask));
grid(loc) = 1:n;
A = sparse(n,n);
for i = 1:size(mask,1)
    for j = 1:size(mask,2)
        if mask(i,j) == 1 & mask_b(i,j) == 0
            ind = grid(i,j);
            A(ind,ind) = -4;
            indup = grid(i-1,j);inddown = grid(i+1,j);indleft = grid(i,j-1);indright = grid(i,j+1);
            A(ind,indup) = 1;A(ind,inddown) = 1;A(ind,indleft) = 1;A(ind,indright) = 1;
        elseif mask(i,j) == 1
            ind=grid(i,j);
            A(ind,ind)=1;
        end
    end
end
C = zeros(n,1);
final = zeros(size(mask,1),size(mask,2),3);
for c=1:3
    S=double(S0(minx:maxx,miny:maxy,c));
    T=double(T0(minx:maxx,miny:maxy,c));
    for i=1:size(S,1)
        for j=1:size(S,2)
            if mask(i,j) == 1 & mask_b(i,j) == 1
                C(grid(i,j)) = T(i,j);
            elseif mask(i,j) == 1 & mask_edge(i,j) == 1
                v1 = -4*S(i,j) + S(i-1,j) +S(i,j-1) + S(i,j+1) + S(i+1,j);
                C(grid(i,j)) = v1;
            elseif mask(i,j) == 1
                C(grid(i,j)) = T(i,j);
            end
        end
    end
    C1=A\C;
    final(:,:,c) = T;
    for i=1:size(S,1)
        for j=1:size(S,2)
            if(grid(i,j) ~= 0)
                final(i,j,c) = C1(grid(i,j));
            end
        end
    end
end

T0(minx:maxx,miny:maxy,:) = uint8(final);
figure,imshow(T0);