function [divR, deltaf, error, interior_rmse] = verify_euler_lagrange_eq(F, V, f)
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
df_on_edges = full([weight.*tmp1(edge_ind) weight.*tmp2(edge_ind)]);
R = zeros(size(E,1),2,2);
for row=1:size(E,1)
    R(row,:,:) = polar_newton( ...
    [real(df_on_edges(row,1)) real(df_on_edges(row,2)); ...
     imag(df_on_edges(row,1)) imag(df_on_edges(row,2))]);
end
R = reshape(complex(R(:,1,:), R(:,2,:)),size(E,1),2);

% compute theta for the closest rotation matrix
%theta = atan((imag(df_on_edges(:,1))-real(df_on_edges(:,2))) ./ ...
%    (real(df_on_edges(:,1))+imag(df_on_edges(:,2))));
% at each vertex R = [cos(theta) -sin(theta); sin(theta) cos(theta)];
%R = [complex(cos(theta), -sin(theta)) complex(sin(theta), cos(theta))];

divR = compute_divergence(F, V, E, cotan, R, true);

deltaf = -L * f;

%rmse_relative = sqrt(mean((abs(divR - deltaf) ./ abs(divR)) .^ 2));
rmse_abs = sqrt(mean((abs(divR - deltaf)) .^ 2));
total_l2 = sqrt(sum((abs(divR - deltaf)) .^ 2));
%disp(['Root mean square relative error: ', num2str(rmse_relative)]);
%disp(['Root mean square absolute error: ', num2str(rmse_abs)]);
disp(['Total l2 error: ', num2str(total_l2)]);
disp(['Total l2 norm of divR: ', num2str(sqrt(sum(abs(divR).^2)))]);
disp(['Total l2 norm of laplacian of f: ', num2str(sqrt(sum(abs(deltaf).^2)))]);

[boundary_vertices, boundary_edge_ind, boundary_normal] = identify_boundary(F, V, E);
interior = setdiff(1:size(V,1), boundary_vertices);
total_l2_interior = sqrt(sum((abs(divR(interior) - deltaf(interior))) .^ 2));
disp(['Total l2 error in the interior: ', num2str(total_l2_interior)]);
disp(['Total l2 norm of divR in the interior: ', num2str(sqrt(sum(abs(divR(interior)).^2)))]);
disp(['Total l2 norm of laplacian of f in the interior: ', num2str(sqrt(sum(abs(deltaf(interior)).^2)))]);

%median_relative = median(abs(divR - deltaf) ./ abs(divR));
%median_abs = median(abs(divR - deltaf));
%disp(['Median relative error: ', num2str(median_relative)]);
%disp(['Median absolute error: ', num2str(median_abs)]);

error = abs(divR - deltaf);
interior_rmse = sqrt(mean(error(interior) .^ 2));

% force the laplacian of f to be 0 for display purpose on the boundary
deltaf(boundary_vertices) = 0;

end