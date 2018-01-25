function [dualh] = dual_map_on_vertices(F, V, primalh)
E = expand_faces_to_edges(F);
cotan = compute_cotan(V, E);

edge_len = twodnorm(V(E(:,1),:) - V(E(:,2),:));
dual_area_contribution = 0.125 * (edge_len .^ 2) .* cotan;
dualarea = accumarray(E(:,1), dual_area_contribution, [size(V,1) 1]);
dualh = zeros(size(V,1), size(primalh,2));
for i=1:size(primalh,2)
    dualh(:,i) = accumarray(E(:,1), ...
        dual_area_contribution .* repmat(primalh(:,i),6,1), [size(V,1) 1]);
end
dualh = dualh ./ dualarea;
end
