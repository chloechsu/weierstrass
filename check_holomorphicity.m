function [ error ] = check_holomorphicity( F, V, f )
% Check the curl of [real(f) -imag(f)] and the curl of [imag(f) real(f)]
E = expand_faces_to_edges(F);
h_on_edges = 0.5 * (f(E(:,1),:) + f(E(:,2),:));
curl1 = compute_curl(F, V, [real(h_on_edges) -imag(h_on_edges)]);
curl2 = compute_curl(F, V, [imag(h_on_edges) real(h_on_edges)]);
disp(['Curl of corresponding vector fields: ', num2str(rmse(curl1)), ...
    ' ', num2str(rmse(curl2)), ...
    ' average: ', num2str(mean([rmse(curl1), rmse(curl2)]))]);
error = 0.5*abs(curl1) + 0.5*abs(curl2);

% Check the relative error in holomorphicity (how far f is from being
% holomorphic)
%[fz, fzbar] = decompose_df(F, V, f);
%rmse_relative = sqrt(mean((abs(fzbar) ./ abs(fz)) .^ 2));
%rmse_abs = sqrt(mean((abs(fzbar)) .^ 2));
%disp(['Root mean square relative error: ', num2str(rmse_relative)]);
%disp(['Root mean square absolute error: ', num2str(rmse_abs)]);
%error = abs(fzbar);
end

function [e] = rmse(error)
e = sqrt(mean((abs(error)) .^ 2));
end