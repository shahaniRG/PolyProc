[gid_map,removedFront,removedEnd] = trim_map(DS.GIDvol,0);

function [gid_map,removedFront,removedEnd] = trim_map(gid_map, airIdx)
% Marcel Chlupsa 
% v0.1(6/22) opened
% v0.3(10/24/22) edited
% v0.5(6/9/25) replaced mean==airIdx with a;;(x==airIdx), replaced updating 
% k with end, added details on how many slices were removed

if ndims(gid_map)~=3
    error('input should be 3-dimensional!\nr')
end

% initialize remove counters
removedFront = [0 0 0];
removedEnd = [0 0 0];

% 3rd dim
while all( gid_map(:,:,end) == airIdx, 'all' )
    gid_map(:,:,end) = [];
    removedEnd(3) = removedEnd(3) + 1;
end
while all( gid_map(:,:,1) == airIdx, 'all' )
    gid_map(:,:,1) = [];
    removedFront(3) = removedFront(3) + 1;
end

% 2nd dim
while all( gid_map(:,end,:) == airIdx, 'all' )
    gid_map(:,end,:) = [];
    removedEnd(2) = removedEnd(2) + 1;
end
while all( gid_map(:,1,:) == airIdx, 'all' )
    gid_map(:,1,:) = [];
    removedFront(2) = removedFront(2) + 1;
end

% 1st dim
while all( gid_map(end,:,:) == airIdx, 'all' )
    gid_map(end,:,:) = [];
    removedEnd(1) = removedEnd(1) + 1;
end
while all( gid_map(1,:,:) == airIdx, 'all' )
    gid_map(1,:,:) = [];
    removedFront(1) = removedFront(1) + 1;
end

% add boundary layer (e.g. for prisms-pf)
a = airIdx*ones(1,size(gid_map,2),size(gid_map,3));
gid_map = cat(1,a,gid_map,a);
b = airIdx*ones(size(gid_map,1),1,size(gid_map,3));
gid_map = cat(2,b,gid_map,b);
c = airIdx*ones(size(gid_map,1),size(gid_map,2),1);
gid_map = cat(3,c,gid_map,c);

% we added a one-layer buffer region so remove that from the removed counts
removedFront = removedFront-1;
removedEnd = removedEnd-1;
end
