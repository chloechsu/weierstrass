function [div] = compute_divergence(F, V, E, cotan, vector_on_edges, closed_boundary)
edgevector = V(E(:,2),:) - V(E(:,1),:);
edge_flow = sum(vector_on_edges .* edgevector, 2);
div = accumarray(E(:,1), 0.5 * cotan .* edge_flow, [size(V,1) 1]);

if(closed_boundary)
    % on the boundary, when computing the flux of dual cells, close the dual
    % cells with primal boundary edges
    [boundary_v, boundary_ind, boundary_normal] = identify_boundary(F, V, E);
    edge_flow_boundary = sum(vector_on_edges(boundary_ind,:) .* boundary_normal ./ 2, 2);
    div = div + accumarray(E(boundary_ind, 1), edge_flow_boundary, [size(V,1) 1]);
end

end

