function [div] = compute_divergence(F, V, E, cotan, vector_on_edges)
edgevector = V(E(:,2),:) - V(E(:,1),:);
edge_flow = sum(vector_on_edges .* edgevector, 2);
div = accumarray(E(:,1), 0.5 * cotan .* edge_flow, [size(V,1) 1]);
end