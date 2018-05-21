% This is a script that demos computing Bounded Biharmonic Weights
% automatically for a 2D shape.
%
% This file and any included files (unless otherwise noted) are copyright Alec
% Jacobson. Email jacobson@inf.ethz.ch if you have questions
%
% Copyright 2011, Alec Jacobson (jacobson@inf.ethz.ch)
%

% NOTE: Please contact Alec Jacobson, jacobson@inf.ethz.ch before
% using this code outside of an informal setting, i.e. for comparisons.

%
results = [];
for i=1:7
    mesh_source = ...
        strcat(['square',num2str(2^i), 'x', num2str(2^i), '.obj']);
    C = [-0.1,-0.1; 1.1,-0.1; 1.1,1.1; -0.1,1.1];
    new_C = [-0.1,-0.1; 1.1,-0.1; 1.1,1.1; -0.5,1.5];
    results = [results; experiment(mesh_source, C, new_C)];
end
csvwrite('square_results_top_left_extend.csv',results); 
%

%
results = [];
for i=1:7
    mesh_source = ...
        strcat(['square',num2str(2^i), 'x', num2str(2^i), '.obj']);
    C = [-0.1,-0.1; 1.1,-0.1; 1.1,1.1; -0.1,1.1];
    new_C = [-0.1,-0.1; 1.1,-0.1; 1.1,1.1; 0.3,0.7];
    results = [results; experiment(mesh_source, C, new_C)];
end
csvwrite('square_results_top_left_shrink.csv',results); 
%

%
results = [];
for i=1:7
    mesh_source = ...
        strcat(['square',num2str(2^i), 'x', num2str(2^i), '.obj']);
    C = [-0.1,-0.1; 1.1,-0.1; 1.1,1.1; -0.1,1.1];
    new_C = [-0.1,-0.1; 1.1,-0.1; 1.1,1.1; 0.1,0.9];
    results = [results; experiment(mesh_source, C, new_C)];
end
csvwrite('square_results_top_left_small_shrink.csv',results); 
%

%
results = [];
for i=1:6
    mesh_source = ...
        strcat(['L-shape_',num2str(2^(i-1)*5), 'x', num2str(2^(i-1)*5), '.obj']);
    C = [-0.05,-0.05; 0.55,-0.05; 1.05,-0.05; 1.05,0.55; ...
         0.55,0.55; 0.55,1.05; -0.05,1.05; -0.05,0.55];
    new_C = [-0.05,-0.05; 0.55,-0.05; 1.05,-0.05; 1.05,0.55; ...
         0.7,0.7; 0.55,1.05; -0.05,1.05; -0.05,0.55];
    results = [results; experiment(mesh_source, C, new_C)];
end
csvwrite('L_results_center_extend.csv',results); 
%

%
results = [];
for i=1:6
    disp(i);
    mesh_source = ...
        strcat(['L-shape_',num2str(2^(i-1)*5), 'x', num2str(2^(i-1)*5), '.obj']);
    C = [-0.05,-0.05; 0.55,-0.05; 1.05,-0.05; 1.05,0.55; ...
         0.55,0.55; 0.55,1.05; -0.05,1.05; -0.05,0.55];
    new_C = [-0.05,-0.05; 0.55,-0.05; 1.05,-0.05; 1.05,0.55; ...
         0.2,0.22; 0.55,1.05; -0.05,1.05; -0.05,0.55];
    results = [results; experiment(mesh_source, C, new_C)];
end
csvwrite('L_results_center_shrink.csv',results); 
%

results = [];
for i=1:6
    disp(i);
    mesh_source = ...
        strcat(['annulus_',num2str(2^i), 'x', num2str(2^i), '.obj']);
    C = [-1.05,-1.05; 1.05,-1.05; 1.05,1.05; -1.05,1.05];
    new_C = [-1.05,-1.05; 1.05,-1.05; -0.5,-0.55; -1.05,1.05];
    results = [results; experiment(mesh_source, C, new_C)];
end
csvwrite('annulus_squish_corner.csv',results); 

function stats_comparison = experiment(mesh_source, C, new_C)

[V,F] = load_mesh(mesh_source);
V = V(:,1:2);
T = new_C - C;
W = cauchy_green_weights(V,C);
 
g = compute_holomorphic_stress(F, V, W, C, T);
[f_weierstrass, branch_points] = compute_elastic_map(F, V, g, false);
[f_polar, branch_points] = compute_elastic_map_polar(F, V, g);

f_cauchygreen = W * complex(T(:,1),T(:,2));

stats_comparison = [get_stats(F, V, g, f_weierstrass), ...
    get_stats(F, V, g, f_polar), get_stats(F, V, g, f_cauchygreen)];
idx = [1 4 7 2 5 8 3 6 9];
stats_comparison = stats_comparison(:,idx);

end

function stats = get_stats(F, V, g, f)
    [fz, fzbar] = decompose_df(F, V, f);
    g_recovered = dual_map_on_vertices(F, V, 2 * fz - fz ./ abs(fz));
    g_diff = sqrt(mean(abs(g - g_recovered) .^ 2));
    [divR, laplacian_of_f, euler_lagrange_error, elas_err] = ...
            verify_euler_lagrange_eq(F, V, f);
    energy = elastic_energy(V, F, fz, fzbar);
    stats = [g_diff, elas_err, energy];
end
