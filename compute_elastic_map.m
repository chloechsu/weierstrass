function [ elastic_f, branch_points ] = compute_elastic_map( F, V, g_on_vertices, conformal )
% Compute the elastic map that minimizes elastic energy from the stress
% tensor g, using the theory of Weierstrass representation
% If conformal = true, returns an elastic map that is also conformal,
% otherwise the returned elastic map may not be conformal (holomorphic).

G = find_antiderivative(F, V, g_on_vertices);

% to avoid sigularity at the branch cut of the complex square root, 
% use the double ranched cover when necessary
[F_, V_, h_, branch_points, is_covering_disconnected] = ...
    double_branch_square_root(F, V, g_on_vertices);
%h_error = check_holomorphicity(F_, V_, h_, 'full');
%if size(F_,1) > size(F,1)
%    h_error = (h_error(1:size(F,1)) + h_error(size(F,1)+1:2*size(F,1))) / 2;
%end

if is_covering_disconnected
    % there are two disconnected sheets, so we need to specify H_ at two
    % vertices
    H_ = find_antiderivative(F_, V_, h_, 1, size(V,1)+1);
elseif size(branch_points,1) == 1
    % shift H by a constant so that H is zero at branch points
    % where h has a zero of degree n, H should have a zero of degree n+1
    H_ = find_antiderivative(F_, V_, h_, size(V,1)*2+1);
    H_ = H_ - H_(2*size(V,1)+1);
elseif size(branch_points,1) == 0
    % the usual case with principal square root without considering branch cut
    H_ = find_antiderivative(F_, V_, h_);
else
    % TODO: what if there are more than 1 branch points?
end

% compute elastic map from h, H, G
% use the "upper layer" of the double branched cover for h and H
H = H_(1:size(V,1));
h = h_(1:size(V,1));

% build k to cancel out the singularity points in H / conj(h)
% assuming all the zeros of h are first order
k = zeros(size(V,1),1);
%for i=1:size(branch_points,1)
%   z0 = complex(branch_points(i,1), branch_points(i,2));
%   H0 = H_(2*size(V,1)+i);
%   k = k - 0.5 * H0 ./ sqrt(complex(V(:,1),V(:,2)) - z0);
%end

elastic_f = 0.5 * (G + H ./ conj(h)) + conj(k);

if(conformal)
    [fz, fz_bar] = decompose_df(F, V, elastic_f);
    fz_on_vertices = dual_map_on_vertices(F, V, fz);
    elastic_f = find_antiderivative(F, V, fz_on_vertices);
end

end