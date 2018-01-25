function [area] = dual_area(F, V)
E = expand_faces_to_edges(F);
cotan = compute_cotan(V, E);

edge_len = twodnorm(V(E(:,1),:) - V(E(:,2),:));
dual_area_contribution = 0.125 * (edge_len .^ 2) .* cotan;
area = accumarray(E(:,1), dual_area_contribution, [size(V,1) 1]);
end
