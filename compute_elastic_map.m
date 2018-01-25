function [ elastic_f ] = compute_elastic_map( F, V, g_on_vertices)
% Compute the elastic map that minimizes elastic energy from the stress
% tensor g, using the theory of Weierstrass representation

% circshift the vertices in each face to have a row for each edge
E = expand_faces_to_edges(F);
cotan = compute_cotan(V, E);
L = compute_laplacian(F, V, E, cotan);

G = find_antiderivative(F, V, E, cotan, L, g_on_vertices);
disp('Verifying holomorphicity for G...')
check_holomorphicity(F, V, G);

h_on_vertices = sqrt(g_on_vertices);
disp('Verifying holomorphicity for h on vertices...')
check_holomorphicity(F, V, h_on_vertices);
H = find_antiderivative(F, V, E, cotan, L, h_on_vertices);
disp('Verifying holomorphicity for H...')
check_holomorphicity(F, V, H);

% compute elastic map from h, H, G
elastic_f = 0.5 * (G + H ./ conj(h_on_vertices));
%disp(size(elastic_f));

%[elastic_fz, elastic_fzbar] = decompose_df(F, V, elastic_f);
%recovered_g = 2 * elastic_fz - elastic_fz ./ abs(elastic_fz);
%disp('Relative error in recovered stress:');
%disp(norm(g-recovered_g) / norm(g));

end

function [G] = find_antiderivative(F, V, E, cotan, L, g_on_vertices)
gradGreal = [real(g_on_vertices) -imag(g_on_vertices)];
gradGimag = [imag(g_on_vertices) real(g_on_vertices)];
divdGreal = compute_divergence(F, V, E, cotan, ...
    0.5*gradGreal(E(:,1),:)+0.5*gradGreal(E(:,2),:));
divdGimag = compute_divergence(F, V, E, cotan, ...
    0.5*gradGimag(E(:,1),:)+0.5*gradGimag(E(:,2),:));
Greal = solve_laplacian(L, divdGreal);
Gimag = solve_laplacian(L, divdGimag);
%disp(sqrt(mean((gradGreal - dual_map_on_vertices(F, V, compute_df(F, V, Greal))) .^ 2)));
%disp(sqrt(mean((gradGimag - dual_map_on_vertices(F, V, compute_df(F, V, Gimag))) .^ 2)));
G = complex(Greal, Gimag);
end

function [x] = solve_laplacian(L, b)
x = zeros([size(L,1) 1]);
x(2:end) = -L(2:end,2:end) \ b(2:end);
end