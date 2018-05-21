function [ elastic_f, branch_points ] = compute_elastic_map_polar( F, V, g_on_vertices )
% Compute the elastic map that minimizes elastic energy from the stress
% alternative approach using g = (trace Y - 1) R and take df = (trace Y/2)R
df = 0.5 * (abs(g_on_vertices)+1) ./ abs(g_on_vertices) .* g_on_vertices;
elastic_f = find_antiderivative(F, V, df);
branch_points = [];
end
