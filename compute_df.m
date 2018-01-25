function [ df ] = compute_df( F, V, f )
% Using discrete gradient, the gradient on each face is (vi, vj, vk) is
% (fj - fi) (vi - vk)^perp / (2 Area) + (fk - fi) (vj - vi)^perp / (2 Area)

vi_vk = V(F(:,1),:) - V(F(:,3),:);
vi_vk_perp = [vi_vk(:,2) -vi_vk(:,1)];
vj_vi = V(F(:,2),:) - V(F(:,1),:);
vj_vi_perp = [vj_vi(:,2) -vj_vi(:,1)];

% signed area
% area is positive when the triangle vi vj vk is oriented clockwise
A = V(F(:,2),1) .* V(F(:,3),2) + V(F(:,3),1) .* V(F(:,1),2) + V(F(:,1),1) .* V(F(:,2),2);
A = A - V(F(:,2),1) .* V(F(:,1),2) - V(F(:,1),1) .* V(F(:,3),2) - V(F(:,3),1) .* V(F(:,2),2);
A = -A / 2;
df = repmat(f(F(:,2)) - f(F(:,1)), [1,2]) .* vi_vk_perp;
df = df + repmat(f(F(:,3)) - f(F(:,1)), [1,2]) .* vj_vi_perp;
df = df ./ repmat(2 * A, [1,2]);
end