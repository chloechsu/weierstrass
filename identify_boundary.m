function [boundary_vertices, boundary_edge_ind, boundary_normal] = identify_boundary(F, V, E)
% return boundary vertices, indices of boundary edges in E, and corresponding outward
% normal vectors for each edge

sz = [size(V,1) size(V,1)];
ind = sub2ind(sz, E(:,1), E(:,2));
triangle_count = accumarray(ind, 1, [prod(sz) 1], [], [], true);
triangle_count = reshape(triangle_count, sz);

E_unflipped = [F; circshift(F,[0,1]); circshift(F,[0,2])];
E_unflipped = E_unflipped(:,1:2);
boundary_edge_ind = find(triangle_count(sub2ind(sz, E_unflipped(:,1), E_unflipped(:,2))) == 1);
% note that the faces are oriented counterclockwise
rotate_90_clockwise = [0 1; -1 0];
edgevector = V(E(boundary_edge_ind,2),:) - V(E(boundary_edge_ind,1),:);
boundary_normal = edgevector * transpose(rotate_90_clockwise);

% include a flipped copy of all edges
boundary_normal = [boundary_normal; boundary_normal];
boundary_edge_ind = [boundary_edge_ind; boundary_edge_ind + size(E_unflipped,1)];
boundary_vertices = unique(E(boundary_edge_ind,1));

end