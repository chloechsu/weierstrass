function [ EV ] = elastic_energy( V, F, fz, fzbar )
% computes the elastic energy of a deformation

A = V(F(:,2),1) .* V(F(:,3),2) + V(F(:,3),1) .* V(F(:,1),2) + V(F(:,1),1) .* V(F(:,2),2);
A = A - V(F(:,2),1) .* V(F(:,1),2) - V(F(:,1),1) .* V(F(:,3),2) - V(F(:,3),1) .* V(F(:,2),2);
A = abs(A) / 2;
EV = 0.5 * sum( ((abs(fz)-1).^2 + abs(fzbar).^2) .* A );

end

