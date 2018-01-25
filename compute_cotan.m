function [cot] = compute_cotan(V, F)
vi_vk = V(F(:,1),:) - V(F(:,3),:);
vj_vk = V(F(:,2),:) - V(F(:,3),:);
cos = dot(vi_vk, vj_vk, 2) ./ twodnorm(vi_vk) ./ twodnorm(vj_vk);
cot = cos ./ sqrt(1 - cos .* cos);
end