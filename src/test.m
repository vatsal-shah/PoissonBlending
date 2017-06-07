R = 'L';
n = 32;
G = numgrid(R,n);
g = numgrid(R,12);
D = -delsq(G);
D = full(D);
grid = g;
mask = g>0;
mask_b = mask - imerode(mask,[0 1 0;1 1 1;0 1 0]);
mask_b=logical(mask_b);

for i = 1:size(mask,1)
    for j = 1:size(mask,2)
        if mask(i,j) == 1 & mask_b(i,j) == 1
            ind = grid(i,j);
            D(ind,:) = 0;
            D(ind,ind) = 1;
        end
    end
end

A = zeros(size(D));
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