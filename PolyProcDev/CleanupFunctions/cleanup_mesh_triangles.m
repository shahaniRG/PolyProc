function [F,V_new] = cleanup_mesh_triangles(F, V, n_iter, alpha, lambda)
    % Smooths mesh vertices to equalize edge lengths while preserving shape.
    %
    % Inputs:
    %   F        - faces (Nx3)
    %   V        - vertices (Mx3)
    %   n_iter   - number of smoothing iterations
    %   alpha    - step size for edge equalization
    %   lambda   - regularization weight for shape preservation
    %
    % Outputs:
    %   V_new    - updated vertex positions
    %   F        - faces (unchanged)

    V_new = V;
    N = size(V, 1);

    % Identify boundary vertices
    TR = triangulation(F, V);
    boundaryEdges = freeBoundary(TR);
    boundaryVerts = unique(boundaryEdges(:));

    % Build vertex adjacency list
    E_all = TR.edges;
    adj_list = cell(N, 1);
    for i = 1:size(E_all, 1)
        adj_list{E_all(i,1)}(end+1) = E_all(i,2);
        adj_list{E_all(i,2)}(end+1) = E_all(i,1);
    end
    for i = 1:N
        adj_list{i} = unique(adj_list{i});
    end

    % Face attachments per vertex
    face_adj = vertexAttachments(TR);  % cell array: {faces attached to vertex i}

    for iter = 1:n_iter
        V_temp = V_new;
        TR = triangulation(F, V_new);  % update for normal computation

        for i = 1:N
            if ismember(i, boundaryVerts)
                continue;  % skip boundary vertices
            end

            nbrs = adj_list{i};
            if numel(nbrs) < 3
                continue;
            end

            % Edge vectors and lengths
            e_vecs = V_new(nbrs,:) - V_new(i,:);
            e_lens = vecnorm(e_vecs, 2, 2);
            l_avg = mean(e_lens);

            % Normalized edge vectors
            e_unit = e_vecs ./ (e_lens + eps);
            delta = sum((l_avg - e_lens) .* e_unit, 1) / numel(nbrs);

            % Estimate vertex normal from neighboring faces
            f_idx = face_adj{i}(:);
            fn = faceNormal(TR, f_idx);
            normal = mean(fn,1);
            normal = normal / (norm(normal)+eps);

            % Project delta to tangent plane
            delta_proj = delta - dot(delta, normal) * normal;

            % Laplacian regularization
            laplace = mean(V_new(nbrs,:), 1) - V_new(i,:);

            % Update vertex
            V_temp(i,:) = V_new(i,:) + alpha * delta_proj + lambda * laplace;
        end
        V_new = V_temp;
    end
end
