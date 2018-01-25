function [error] = verify_euler_lagrange_eq(F, V, f)
E = expand_faces_to_edges(F);
cotan = compute_cotan(V, E);
L = compute_laplacian(F, V, E, cotan);

df_on_faces = compute_df(F, V, f);
df_on_faces = repmat(df_on_faces,6,1);
edge_ind = sub2ind([size(V,1),size(V,1)],E(:,1),E(:,2));
% the weight is 0.5 for internal edges and 1.0 for boundary edges
weight = accumarray([E(:,1) E(:,2)],1,[],[],[],true);
weight = ones(size(E,1),1) ./ weight(edge_ind);
tmp1 = accumarray([E(:,1) E(:,2)],df_on_faces(:,1),[],[],[],true);
tmp2 = accumarray([E(:,1) E(:,2)],df_on_faces(:,2),[],[],[],true);
df_on_edges = [weight.*tmp1(edge_ind) weight.*tmp2(edge_ind)];
R = zeros(size(E,1),2,2);
for row=1:size(E,1)
R(row,:,:) = polar_newton( ...
    [real(df_on_edges(row,1)) real(df_on_edges(row,2)); ...
     imag(df_on_edges(row,1)) imag(df_on_edges(row,2))]);
end
%df_matrix = zeros([size(E,1) 2 2]);
%df_matrix(:,1,1) = real(df_on_edges(:,1));
%df_matrix(:,1,2) = real(df_on_edges(:,2));
%df_matrix(:,2,1) = imag(df_on_edges(:,1));
%df_matrix(:,2,2) = imag(df_on_edges(:,2));
%df_matrix_cell = mat2cell(df_matrix, ones(1,size(E,1)), [2], [2]);
%[R, P, k] = cellfun(@(A) polar_newton(reshape(A,2,2)), df_matrix_cell, ...
%    'UniformOutput', false);
%R = permute(cat(3, R{:}), [3,1,2]);
%R = reshape(complex(R(:,1,:), R(:,2,:)),size(E,1),2);

R = reshape(complex(R(:,1,:), R(:,2,:)),size(E,1),2);
divR = compute_divergence(F, V, E, cotan, R);

deltaf = -L * f;

rmse_relative = sqrt(mean((abs(divR - deltaf) ./ abs(divR)) .^ 2));
rmse_abs = sqrt(mean((abs(divR - deltaf)) .^ 2));
disp(['Root mean square relative error: ', num2str(rmse_relative)]);
disp(['Root mean square absolute error: ', num2str(rmse_abs)]);

median_relative = median(abs(divR - deltaf) ./ abs(divR));
median_abs = median(abs(divR - deltaf));
disp(['Median relative error: ', num2str(median_relative)]);
disp(['Median absolute error: ', num2str(median_abs)]);

error = abs(divR - deltaf) ./ abs(divR);

disp(mean(abs(divR)));

end