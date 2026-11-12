%% Run Script
cs = crystalSymmetry('cubic');

t1 = import_from_h5("DATA\t1.h5");
t2 = import_from_h5("DATA\t2.h5");

t1_fp = compute_feature_properties(t1);
t2_fp = compute_feature_properties(t2);

t1_fpgn = compute_grain_neighbors(t1_fp);
t2_fpgn = compute_grain_neighbors(t2_fp);

t1_fpgn_cleancomp = cleanup_by_completeness(t1_fpgn);
t2_fpgn_cleancomp = cleanup_by_completeness(t2_fpgn);
t1_fpgn_cleanori = cleanup_by_orientation(t1_fpgn_cleancomp);
t2_fpgn_cleanori = cleanup_by_orientation(t2_fpgn_cleancomp);
t1_fpgn_clean = cleanup_by_volume(t1_fpgn_cleanori);
t2_fpgn_clean = cleanup_by_volume(t2_fpgn_cleanori);

matches = track_grains_mtex( ...
    t1_fpgn_clean.featureData.id, t1_fpgn_clean.featureData.numElement, t1_fpgn_clean.featureData.adjacency, t1_fpgn_clean.featureData.grain_rodV, ...
    t2_fpgn_clean.featureData.id, t2_fpgn_clean.featureData.numElement, t2_fpgn_clean.featureData.adjacency, t2_fpgn_clean.featureData.grain_rodV, cs);

%% Debugging
report = check_matching(matches, t1_fpgn_clean.featureData.numElement(:,1), t2_fpgn_clean.featureData.numElement(:,1));

%% ======================= Grain tracking =======================
function matches = track_grains_mtex(id1, numElement1, adj1, rodV1, ...
                                     id2, numElement2, adj2, rodV2, cs, varargin)
% TRACK_GRAINS_MTEX  Symmetry-aware grain tracking between two states.
%   1) Bilateral top-K orientation seeding among the largest grains
%   2) Neighbor-consensus expansion (one-to-one core)
%   3) Optional min-cost assignment of remaining free grains
%   4) Fragment attachment: a free grain joins a matched grain if it is
%      adjacent to that grain's partner (split) or to a grain already mapped
%      to the same partner (merge), passes neighbor consensus, and is within
%      AttachMaxDeg. Repeats until no change.
%
% Output:
%   matches : [K,2] (state-1 ID, state-2 ID). IDs may repeat in either
%             column (1:N splits, N:1 merges). Use check_matching to group.

p = inputParser;
p.addParameter('SeedThresholdDeg',   2.0,  @(x)isnumeric(x)&&isscalar(x)&&x>0);
p.addParameter('TopNSeeds',          300,  @(x)isnumeric(x)&&isscalar(x)&&x>0);
p.addParameter('TopK',               3,    @(x)isnumeric(x)&&isscalar(x)&&x>=1);
p.addParameter('ExpandThresholdDeg', 6.0,  @(x)isnumeric(x)&&isscalar(x)&&x>0);
p.addParameter('MinNbrVotes',        1,    @(x)isnumeric(x)&&isscalar(x)&&x>=0);
p.addParameter('MaxRounds',          12,   @(x)isnumeric(x)&&isscalar(x)&&x>=1);
p.addParameter('DoAssignment',       true, @(x)islogical(x)&&isscalar(x));
p.addParameter('AssignMaxDeg',       8.0,  @(x)isnumeric(x)&&isscalar(x)&&x>0);
p.addParameter('SizeWeight',         0.1,  @(x)isnumeric(x)&&isscalar(x)&&x>=0);
p.addParameter('AttachMaxDeg',       2.0,  @(x)isnumeric(x)&&isscalar(x)&&x>0);
p.parse(varargin{:});
par = p.Results;

% ---------- inputs (all work below is in index space) ----------
ids1 = double(id1(:,1));  sz1 = double(numElement1(:,1));
ids2 = double(id2(:,1));  sz2 = double(numElement2(:,1));
N1 = numel(ids1);  N2 = numel(ids2);

rodV1 = ensure3xN(rodV1);  rodV2 = ensure3xN(rodV2);
assert(size(rodV1,2) == N1, 'rodV1 must have one column per row of numElement1.');
assert(size(rodV2,2) == N2, 'rodV2 must have one column per row of numElement2.');
q1 = rod2quat(rodV1.');
q2 = rod2quat(rodV2.');
o1 = orientation(rotation(q1(1),q1(2),q1(3),q1(4)), cs);
o2 = orientation(rotation(q2(1),q2(2),q2(3),q2(4)), cs);

A1 = adj_sparse(adj1, ids1);  nb1 = adj_lists(A1);   % N1 x N1, symmetric
A2 = adj_sparse(adj2, ids2);  nb2 = adj_lists(A2);   % N2 x N2, symmetric

% ---------- 1) bilateral top-K seeding ----------
s1 = top_n(sz1, par.TopNSeeds);
s2 = top_n(sz2, par.TopNSeeds);
D12 = angle_outer(o1(s1), o2) / degree;              % |s1| x N2
D21 = angle_outer(o2(s2), o1) / degree;              % |s2| x N1
K12 = topk_mask(D12, par.TopK, par.SeedThresholdDeg);
K21 = topk_mask(D21, par.TopK, par.SeedThresholdDeg);
[a, b] = find(K12(:, s2) & K21(:, s1).');            % mutual top-K
si = s1(a);  sj = s2(b);
sang = D12(sub2ind(size(D12), a(:), sj(:)));

p12 = zeros(N1,1);  p21 = zeros(N2,1);               % 1:1 core partners (0 = free)
[~, ord] = sort(sang);                               % lowest misorientation first
for t = ord(:).'
    if p12(si(t)) == 0 && p21(sj(t)) == 0
        p12(si(t)) = sj(t);  p21(sj(t)) = si(t);
    end
end
if ~any(p12), matches = zeros(0,2); return; end

% ---------- 2) neighbor-consensus expansion ----------
for rnd = 1:par.MaxRounds
    nAdded = 0;
    for i = find(p12).'
        j  = p12(i);
        c1 = nb1{i}(p12(nb1{i}) == 0);
        c2 = nb2{j}(p21(nb2{j}) == 0);
        if isempty(c1) || isempty(c2), continue; end
        D = angle_outer(o1(c1), o2(c2)) / degree;
        [dv, lin] = sort(D(:));
        for t = 1:numel(lin)
            if dv(t) > par.ExpandThresholdDeg, break; end
            [a, b] = ind2sub(size(D), lin(t));
            ii = c1(a);  jj = c2(b);
            if p12(ii) || p21(jj), continue; end
            if nbr_votes(p12, nb1, nb2, ii, jj) >= par.MinNbrVotes
                p12(ii) = jj;  p21(jj) = ii;  nAdded = nAdded + 1;
            end
        end
    end
    if nAdded == 0, break; end
end

% ---------- 3) min-cost assignment of free grains ----------
if par.DoAssignment
    f1 = find(p12 == 0);  f2 = find(p21 == 0);
    if ~isempty(f1) && ~isempty(f2)
        M = sparse(find(p12), p12(p12 > 0), 1, N1, N2);
        V = A1(f1,:) * M * A2(:,f2);                 % matched neighbor pairs
        D = angle_outer(o1(f1), o2(f2)) / degree;
        [a, b] = find(D <= par.AssignMaxDeg & V >= par.MinNbrVotes);
        if ~isempty(a)
            ii = f1(a(:));  jj = f2(b(:));
            c  = D(sub2ind(size(D), a(:), b(:))) + ...
                 par.SizeWeight * abs(log((sz2(jj) + eps) ./ (sz1(ii) + eps)));
            [r, q] = sparse_assign(ii, jj, c, N1, N2);
            p12(r) = q;  p21(q) = r;
        end
    end
end

% ---------- 4) fragment attachment (1:N splits, N:1 merges) ----------
P = [find(p12), p12(p12 > 0)];
for rnd = 1:par.MaxRounds
    M  = sparse(P(:,1), P(:,2), 1, N1, N2);
    u1 = find(~any(M, 2));
    u2 = find(~any(M, 1)).';
    add = zeros(0,2);

    % split: free state-2 grain j is adjacent to a partner of matched grain i
    if ~isempty(u2)
        S = M * A2(:, u2);                           % # partners of i adjacent to j
        F = (S > 0) & (A1 * S >= par.MinNbrVotes);   % + neighbor consensus
        for k = find(any(F, 1))
            j = u2(k);  c = find(F(:, k));
            [d, t] = min(angle(o1(c), o2(j)) / degree);
            if d <= par.AttachMaxDeg, add(end+1,:) = [c(t), j]; end %#ok<AGROW>
        end
    end

    % merge: free state-1 grain i is adjacent to a grain already mapped to j
    if ~isempty(u1)
        S = A1(u1, :) * M;                           % # neighbors of i mapped to j
        F = (S > 0) & (S * A2 >= par.MinNbrVotes);
        for k = find(any(F, 2)).'
            i = u1(k);  c = find(F(k, :)).';
            [d, t] = min(angle(o2(c), o1(i)) / degree);
            if d <= par.AttachMaxDeg, add(end+1,:) = [i, c(t)]; end %#ok<AGROW>
        end
    end

    if isempty(add), break; end
    P = [P; add]; %#ok<AGROW>
end

% ---------- output, ordered by size ----------
[~, ord] = sortrows([-sz1(P(:,1)), -sz2(P(:,2))]);
P = P(ord, :);
matches = [ids1(P(:,1)), ids2(P(:,2))];
end

% ======================= grain helpers =======================
function rod = ensure3xN(rod)
    rod = double(rod);
    if size(rod,1) == 3
        % ok
    elseif size(rod,2) == 3
        rod = rod.';
    else
        error('rodV must be 3xN or Nx3');
    end
end

function A = adj_sparse(adj, ids)
    % Symmetric index-space adjacency; drops self-loops and unknown IDs.
    n = numel(ids);
    if isempty(adj), A = sparse(n, n); return; end
    [ok1, a] = ismember(double(adj(:,1)), ids);
    [ok2, b] = ismember(double(adj(:,2)), ids);
    k = ok1 & ok2 & a ~= b;
    A = double(sparse([a(k); b(k)], [b(k); a(k)], 1, n, n) > 0);
end

function nb = adj_lists(A)
    n = size(A, 2);
    if nnz(A) == 0, nb = repmat({zeros(0,1)}, n, 1); return; end
    [r, c] = find(A);
    nb = accumarray(c, r, [n 1], @(x){sort(x)});
end

function idx = top_n(sizes, N)
    [~, ord] = sort(sizes, 'descend');
    idx = ord(1:min(N, numel(ord)));
end

function K = topk_mask(D, k, thr)
    % Logical sparse mask of the k smallest entries per row that are <= thr.
    [v, j] = mink(D, k, 2);
    r  = repmat((1:size(D,1)).', 1, size(v,2));
    ok = v <= thr;
    K  = sparse(r(ok), j(ok), true, size(D,1), size(D,2));
end

function v = nbr_votes(p12, nb1, nb2, i, j)
    % # neighbors of i whose 1:1 partner is a neighbor of j
    v = nnz(ismember(p12(nb1{i}), nb2{j}));
end

% ======================= shared =======================
function [r, c] = sparse_assign(ci, cj, cost, n1, n2)
% Min-cost one-to-one assignment over sparse candidates, solved per connected
% component with matchpairs. Infeasible pairs are never returned.
    r = zeros(0,1);  c = zeros(0,1);
    if isempty(ci), return; end
    ci = ci(:);  cj = cj(:);  cost = cost(:);
    comp = conncomp(graph(ci, n1 + cj, ones(numel(ci),1), n1 + n2));
    [~, ~, g] = unique(comp(ci));
    groups = accumarray(g(:), (1:numel(ci)).', [], @(x){x});
    for k = 1:numel(groups)
        e = groups{k};
        [ur, ~, li] = unique(ci(e));
        [uc, ~, lj] = unique(cj(e));
        cU   = max(cost(e)) + 1;                     % cost of leaving a node free
        bigC = 2*cU + 1;                             % worse than leaving both free
        Ck = repmat(bigC, numel(ur), numel(uc));
        Ck(sub2ind(size(Ck), li, lj)) = cost(e);
        pr = matchpairs(Ck, cU);
        keep = Ck(sub2ind(size(Ck), pr(:,1), pr(:,2))) < bigC;
        r = [r; ur(pr(keep,1))]; %#ok<AGROW>
        c = [c; uc(pr(keep,2))]; %#ok<AGROW>
    end
end

%% ======================= Particle tracking =======================
function pmatches = track_particles_by_grain_context( ...
        grain_matches, pgAdj1, pgAdj2, partIDs1, C1, partIDs2, C2, varargin)
% TRACK_PARTICLES_BY_GRAIN_CONTEXT
%   1) Candidates: particle pairs sharing >= 1 mapped grain neighbor
%   2) Global translation (or rigid) fit from high-Jaccard candidates
%   3) One-to-one core by min-cost assignment (distance + Jaccard)
%   4) Attachment of free particles to matched ones using containment
%      overlap |A and B| / min(|A|,|B|), so a fragment whose grain neighbors
%      are a subset of its parent's still qualifies -> 1:N and N:1
%
% grain_matches may be many-to-many (output of track_grains_mtex).
% Grain IDs must be positive integers.
%
% Name-value parameters:
%   'MinNbrCount'     (1)     min # mapped grain neighbors for a state-1 particle
%   'MinOverlap'      (0.5)   Jaccard threshold (core) / containment threshold (attach)
%   'MaxDist'         (12.0)  max distance after transform (centroid units)
%   'TransEstOverlap' (0.8)   Jaccard needed to use a pair in the transform fit
%   'UseRigid'        (false) weighted Kabsch (R,t) instead of translation only
%   'WeightOverlap'   (1.0)   cost = dist/MaxDist + WeightOverlap*(1 - overlap)
%
% Output:
%   pmatches : [K,2] (state-1 particle ID, state-2 particle ID)

p = inputParser;
p.addParameter('MinNbrCount',     1,     @(x)isnumeric(x)&&isscalar(x)&&x>=0);
p.addParameter('MinOverlap',      0.5,   @(x)isnumeric(x)&&isscalar(x)&&x>=0&&x<=1);
p.addParameter('MaxDist',         12.0,  @(x)isnumeric(x)&&isscalar(x)&&x>0);
p.addParameter('TransEstOverlap', 0.8,   @(x)isnumeric(x)&&isscalar(x)&&x>=0&&x<=1);
p.addParameter('UseRigid',        false, @(x)islogical(x)&&isscalar(x));
p.addParameter('WeightOverlap',   1.0,   @(x)isnumeric(x)&&isscalar(x)&&x>=0);
p.parse(varargin{:});
par = p.Results;

partIDs1 = double(partIDs1(:));  partIDs2 = double(partIDs2(:));
C1 = double(C1);  C2 = double(C2);
n1 = numel(partIDs1);  n2 = numel(partIDs2);
assert(size(C1,1) == n1 && size(C2,1) == n2, 'Centroid rows must match particle IDs.');
pmatches = zeros(0,2);
if isempty(pgAdj1) || isempty(pgAdj2) || isempty(grain_matches), return; end

% ---------- 1) grain-context overlap ----------
gm = double(grain_matches);
g1max = max([double(pgAdj1(:,2)); gm(:,1)]);
g2max = max([double(pgAdj2(:,2)); gm(:,2)]);
B1 = incidence(pgAdj1, partIDs1, g1max);             % n1 x g1max
B2 = incidence(pgAdj2, partIDs2, g2max);             % n2 x g2max
Gm = sparse(gm(:,1), gm(:,2), 1, g1max, g2max);      % N:M grain map
MB = double(B1 * Gm > 0);                            % mapped grain neighbors of p1
s1 = full(sum(MB, 2));  s2 = full(sum(B2, 2));
[ci, cj, nI] = find(MB * B2.');                      % shared mapped grains
keep = s1(ci) >= par.MinNbrCount;
ci = ci(keep);  cj = cj(keep);  nI = nI(keep);
if isempty(ci), return; end
jac  = nI ./ (s1(ci) + s2(cj) - nI);
cont = nI ./ min(s1(ci), s2(cj));

% ---------- 2) global translation / rigid fit ----------
[~, o] = sort(jac, 'descend');
[~, f] = unique(ci(o), 'stable');                    % best candidate per p1
b = o(f);
b = b(jac(b) >= par.TransEstOverlap);
R = eye(3);  t = zeros(1,3);
if ~isempty(b)
    Pc = C1(ci(b),:);  Qc = C2(cj(b),:);  W = jac(b);
    if par.UseRigid && numel(b) >= 3
        [R, t] = kabsch_weighted(Pc, Qc, W);
    else
        dd = Qc - Pc;
        t = [wmed(dd(:,1),W), wmed(dd(:,2),W), wmed(dd(:,3),W)];
    end
end
d = vecnorm((C1(ci,:)*R.' + t) - C2(cj,:), 2, 2);

% ---------- 3) one-to-one core ----------
ok   = d <= par.MaxDist & jac >= par.MinOverlap;
cost = d/par.MaxDist + par.WeightOverlap*(1 - jac);
[a1, a2] = sparse_assign(ci(ok), cj(ok), cost(ok), n1, n2);
P = [a1, a2];

% ---------- 4) attachment (1:N splits, N:1 merges) ----------
m1 = false(n1,1);  m1(a1) = true;
m2 = false(n2,1);  m2(a2) = true;
ok   = d <= par.MaxDist & cont >= par.MinOverlap;
cost = d/par.MaxDist + par.WeightOverlap*(1 - cont);

k = find(ok & m1(ci) & ~m2(cj));                     % split: free p2 -> matched p1
[jj, ii] = best_per_key(cj(k), ci(k), cost(k));
P = [P; ii, jj];

k = find(ok & ~m1(ci) & m2(cj));                     % merge: free p1 -> matched p2
[ii, jj] = best_per_key(ci(k), cj(k), cost(k));
P = [P; ii, jj];

pmatches = [partIDs1(P(:,1)), partIDs2(P(:,2))];
end

% ======================= particle helpers =======================
function B = incidence(pgAdj, pids, gmax)
    [ok, r] = ismember(double(pgAdj(:,1)), pids);
    g  = double(pgAdj(:,2));
    ok = ok & g >= 1 & g <= gmax;
    B  = double(sparse(r(ok), g(ok), 1, numel(pids), gmax) > 0);
end

function [key, val] = best_per_key(key, val, c)
    % For each unique key, keep the val with minimum cost.
    [~, o] = sort(c);
    key = key(o);  val = val(o);
    [key, f] = unique(key, 'stable');
    val = val(f);
end

function m = wmed(x, w)
    if any(w <= 0), m = median(x); return; end
    [xs, idx] = sort(x);
    c = cumsum(w(idx)) / sum(w);
    m = xs(find(c >= 0.5, 1, 'first'));
end

function [R, t] = kabsch_weighted(P, Q, w)
    % Q ~ P*R.' + t
    w = w(:) / max(eps, sum(w));
    cP = sum(P.*w, 1);  cQ = sum(Q.*w, 1);
    H = ((P - cP).' .* w.') * (Q - cQ);
    [U, ~, V] = svd(H);
    R = V*U.';
    if det(R) < 0
        V(:,3) = -V(:,3);
        R = V*U.';
    end
    t = cQ - cP*R.';
end

%% ======================= Check =======================
function report = check_matching(matches, ids1, ids2)
% CHECK_MATCHING  Validate IDs and classify match groups.
%   A group is a connected component of the bipartite match graph:
%   1:1, 1:N (split), N:1 (merge), N:M.
ids1 = double(ids1(:));  ids2 = double(ids2(:));
report = struct('num_pairs', size(matches,1));
if isempty(matches), warning('matches is empty.'); return; end
L = double(matches(:,1));  R = double(matches(:,2));

report.ids1_not_in_universe = unique(L(~ismember(L, ids1)));
report.ids2_not_in_universe = unique(R(~ismember(R, ids2)));
[~, first] = unique([L R], 'rows', 'stable');
report.duplicate_rows = setdiff((1:numel(L)).', first);
report.coverage_left  = numel(unique(L)) / max(numel(ids1), 1);
report.coverage_right = numel(unique(R)) / max(numel(ids2), 1);

[uL, ~, a] = unique(L);  [uR, ~, b] = unique(R);
nL = numel(uL);
comp = conncomp(graph(a, nL + b, ones(numel(a),1), nL + numel(uR)));
cL = comp(1:nL).';  cR = comp(nL+1:end).';
nG = max(comp);
nl = accumarray(cL, 1, [nG 1]);  nr = accumarray(cR, 1, [nG 1]);
type = repmat("N:M", nG, 1);
type(nl == 1 & nr == 1) = "1:1";
type(nl == 1 & nr > 1)  = "1:N";
type(nl > 1  & nr == 1) = "N:1";
report.groups = table((1:nG).', type, ...
    accumarray(cL, uL, [nG 1], @(x){sort(x).'}), ...
    accumarray(cR, uR, [nG 1], @(x){sort(x).'}), ...
    'VariableNames', {'group', 'type', 'ids1', 'ids2'});
report.row_group = comp(a).';                        % group of each row of matches
report.ok = isempty(report.ids1_not_in_universe) && ...
            isempty(report.ids2_not_in_universe) && isempty(report.duplicate_rows);

fprintf('\n=== Matching check ===\n');
fprintf('Pairs: %d | Groups: %d\n', report.num_pairs, nG);
fprintf('Coverage: left %.1f%% | right %.1f%%\n', ...
    100*report.coverage_left, 100*report.coverage_right);
fprintf('1:1 %d | 1:N %d | N:1 %d | N:M %d\n', nnz(type == "1:1"), ...
    nnz(type == "1:N"), nnz(type == "N:1"), nnz(type == "N:M"));
if ~isempty(report.ids1_not_in_universe)
    fprintf('WARNING: %d left IDs not in ids1\n', numel(report.ids1_not_in_universe));
end
if ~isempty(report.ids2_not_in_universe)
    fprintf('WARNING: %d right IDs not in ids2\n', numel(report.ids2_not_in_universe));
end
if ~isempty(report.duplicate_rows)
    fprintf('WARNING: %d duplicate rows\n', numel(report.duplicate_rows));
end
if report.ok, fprintf('OK\n'); else, fprintf('NOT OK\n'); end
fprintf('======================\n\n');
end