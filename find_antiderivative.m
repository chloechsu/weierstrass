function [G] = find_antiderivative(F, V, g_on_vertices, varargin)
% varargin gives the points where the antiderivative should be zero
% if it's empty, we force the antiderivative to be zero at the first vertex
if isempty(varargin)
    G_zeros = [1];
else
    G_zeros = cell2mat(varargin)';
end

% circshift the vertices in each face to have a row for each edge
E = expand_faces_to_edges(F);
cotan = compute_cotan(V, E);
L = compute_laplacian(F, V, E, cotan);

gradGreal = [real(g_on_vertices) -imag(g_on_vertices)];
gradGimag = [imag(g_on_vertices) real(g_on_vertices)];
divdGreal = compute_divergence(F, V, E, cotan, ...
    0.5*gradGreal(E(:,1),:)+0.5*gradGreal(E(:,2),:), false);
divdGimag = compute_divergence(F, V, E, cotan, ...
    0.5*gradGimag(E(:,1),:)+0.5*gradGimag(E(:,2),:), false);

Greal = solve_laplacian(L, divdGreal, G_zeros);
Gimag = solve_laplacian(L, divdGimag, G_zeros);

%disp(sqrt(mean((gradGreal - dual_map_on_vertices(F, V, compute_df(F, V, Greal))) .^ 2)));
%disp(sqrt(mean((gradGimag - dual_map_on_vertices(F, V, compute_df(F, V, Gimag))) .^ 2)));
G = complex(Greal, Gimag);
end