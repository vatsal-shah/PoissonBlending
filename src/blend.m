function [ T0 ] = blend( T0,S0,msk )
%BLEND Seamless blending using poisson image blending
% S0 : Image containing selected region of source pasted into desired area
% of target, rest everything black
% T0 : Target image
% msk : Mask

mask0 = logical(msk);

[x,y] = find(mask0);
minx = min(x) - 1;
miny = min(y) - 1;
maxx = max(x) + 1;
maxy = max(y) + 1;
n = length(x);

mask = mask0(minx:maxx,miny:maxy);
mask_b = mask - imerode(mask,[0 1 0;1 1 1;0 1 0]);
mask_b=logical(mask_b);
loc = find(mask(:));
grid = zeros(size(mask));
grid(loc) = 1:n;
%A = sparse(n,n);
cnt = 1;
for i = 1:size(mask,1)
    for j = 1:size(mask,2)
        if mask(i,j) == 1 & mask_b(i,j) == 0
            ind = grid(i,j);
            %A(ind,ind) = -4;
            indup = grid(i-1,j);inddown = grid(i+1,j);indleft = grid(i,j-1);indright = grid(i,j+1);
            %A(ind,indup) = 1;A(ind,inddown) = 1;A(ind,indleft) = 1;A(ind,indright) = 1;
            row(cnt:cnt+4) = ind;
            col(cnt) = ind;col(cnt+1) = indup;col(cnt+2) = inddown;col(cnt+3) = indleft;col(cnt+4) = indright;
            val(cnt) = -4;val(cnt+1:cnt+4) = 1;
            cnt = cnt + 5;
        elseif mask(i,j) == 1
            ind=grid(i,j);
            %A(ind,ind)=1;
            row(cnt) = ind;col(cnt) = ind;val(cnt) = 1;
            cnt = cnt + 1;
        end
    end
end
A = sparse(row,col,val);
%pause;
C = zeros(n,1);
final = zeros(size(mask,1),size(mask,2),3);
for c=1:3
    S=double(S0(minx:maxx,miny:maxy,c));
    T=double(T0(minx:maxx,miny:maxy,c));
    for i=1:size(S,1)
        for j=1:size(S,2)
            if mask(i,j) ==1 & mask_b(i,j)==0
                v1 = -4*S(i,j) + S(i-1,j) +S(i,j-1) + S(i,j+1) + S(i+1,j);
                C(grid(i,j))=v1;
            elseif mask(i,j)==1 
                C(grid(i,j))=T(i,j);
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

end