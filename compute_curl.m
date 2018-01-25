function [curl] = compute_curl(F, V, vector_on_edges)
E = expand_faces_to_edges(F);
edges = V(E(:,2),:) - V(E(:,1),:);
edge_flow = sum(vector_on_edges .* edges, 2);
nf = size(F,1);
circumference = twodnorm(edges(1:nf,:)) + twodnorm(edges(nf+1:2*nf,:)) + ...
    twodnorm(edges(2*nf+1:3*nf,:));
curl = edge_flow(1:nf,:) + edge_flow(nf+1:2*nf,:) + edge_flow(2*nf+1:3*nf,:);
curl = curl ./ circumference;
end