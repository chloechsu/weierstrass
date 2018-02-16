function [ g ] = compute_holomorphic_stress( F, V, W, C, T )
% compute a holormophic function g to represent stress

center = mean(C, 1);
directions = C - repmat(center,[size(C,1),1]);

scaling = 1 + sum(T .* directions, 2) ./ (twodnorm(directions) .^ 2);

new_directions = directions + T;
cross_prod = directions(:,1).*new_directions(:,2) - ...
    new_directions(:,1).*directions(:,2);
rotation_sin = cross_prod ./ twodnorm(directions) ./ twodnorm(new_directions);
rotation_sin(isnan(rotation_sin))=0;
rotation_cos = sqrt(1 - rotation_sin .* rotation_sin);

g_on_controls = complex((2*scaling-1) .* rotation_cos, (2*scaling-1) .* rotation_sin);
g = W * g_on_controls;

end

