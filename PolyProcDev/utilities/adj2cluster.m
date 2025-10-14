function C = adj2cluster(A)

% Symmetrize adjacency matrix
S = A + A';

% Reverse Cuthill-McKee ordering
r = fliplr(symrcm(S));

% Get the clusters
C = {r(1)};
for i = 2:numel(r)
    if any(S(C{end}, r(i)))
        C{end}(end+1) = r(i);
    else
        C{end+1} = r(i);
    end
end